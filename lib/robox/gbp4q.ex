defmodule Robox.Gbp4q do
  @moduledoc """
  Select the most preferable bin from a list of bins that can fit all or most
  items, mirroring `gbp4d_solver_dpp_filt` in gbp/src/gbp4d.cpp.

  The bin list should be sorted by volume so that the first entry is the most
  preferred smallest bin and the last entry the least preferred largest bin.

  Result struct:

    * `p`  - profit vector used for the selected bin
    * `it` - one entry per item, 8-tuple `{x, y, z, w, l, d, h, wt}`
    * `bn` - the bin list, one `{l, d, h, w}` tuple per bin
    * `k`  - selection indicator per item, 0 | 1
    * `f`  - selection indicator per bin, exactly one 1
    * `o`  - objective: total volume of fitted items
    * `ok` - did all items fit into the selected bin?
  """

  alias Robox.Gbp4d

  defstruct [:p, :it, :bn, :k, :f, :o, :ok]

  @type t :: %__MODULE__{}

  @tol 1.0e-8

  @doc """
  Pack items `ldhw` (list of `[l, d, h, w]` per item) selecting the most
  preferable bin from `bns` (list of `[l, d, h, w]` per bin).
  """
  @spec pack([[number]], [[number]]) :: t
  def pack(ldhw, bns) do
    ldhw = Enum.map(ldhw, fn [l, d, h, w] -> {l / 1, d / 1, h / 1, w / 1} end)
    bns = Enum.map(bns, fn [l, d, h, w] -> {l / 1, d / 1, h / 1, w / 1} end)

    gbp4d_solver_dpp_filt(ldhw, bns)
  end

  @doc """
  Solve gbp4d w.r.t selecting the most preferable, often smallest, bin from
  the bin list. `ldhw` is a list of `{l, d, h, w}` item tuples and `bns` a
  list of `{l, d, h, w}` bin tuples.
  """
  @spec gbp4d_solver_dpp_filt([tuple], [tuple]) :: t
  def gbp4d_solver_dpp_filt(ldhw, bns) when length(bns) > 0 do
    nbn = length(bns)

    # per item (l, d, h) sorted descending, weight kept aside
    sldhw = Enum.map(ldhw, fn {l, d, h, w} -> {Enum.sort([l, d, h], :desc), w} end)

    # total volume and weight of all items
    vldhw = {
      Enum.reduce(ldhw, 0.0, fn {l, d, h, _}, acc -> acc + l * d * h end),
      Enum.reduce(ldhw, 0.0, fn {_, _, _, w}, acc -> acc + w end)
    }

    sm = Enum.map(bns, fn {l, d, h, w} -> {Enum.sort([l, d, h], :desc), w} end) |> List.to_tuple()
    vm = Enum.map(bns, fn {l, d, h, w} -> {l * d * h, w} end) |> List.to_tuple()

    bnt = List.to_tuple(bns)

    flmt = Enum.to_list(0..(nbn - 1))

    # corner case: a single item fast filt
    flmt =
      if length(ldhw) == 1 do
        filt_fast(sldhw, vldhw, sm, vm, flmt)
      else
        flmt
      end

    # init fit into the last (largest) bin candidate to drive down search space
    {id, flmt} =
      case flmt do
        [] -> {nbn - 1, []}
        _ -> {List.last(flmt), Enum.drop(flmt, -1)}
      end

    p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, elem(bnt, id))
    fit = Gbp4d.gbp4d_solver_dpp(p, ldhw, elem(bnt, id))

    flmt =
      if fit.ok do
        filt_fast(sldhw, vldhw, sm, vm, flmt)
      else
        filt_slow(fit.o, vm, flmt)
      end

    # fit into each remaining bin candidate via approx binary search
    {id, _p, fit} = search(ldhw, bnt, sldhw, vldhw, sm, vm, flmt, id, p, fit)

    f = Enum.map(0..(nbn - 1), fn i -> if i == id, do: 1, else: 0 end)

    %__MODULE__{p: fit.p, it: fit.it, bn: bns, k: fit.k, f: f, o: fit.o, ok: fit.ok}
  end

  defp search(_ldhw, _bnt, _sldhw, _vldhw, _sm, _vm, [], id, p, fit), do: {id, p, fit}

  defp search(ldhw, bnt, sldhw, vldhw, sm, vm, flmt, id, p, fit) do
    id0 = Enum.at(flmt, div(length(flmt), 2))

    p0 = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, elem(bnt, id0))
    fit0 = Gbp4d.gbp4d_solver_dpp(p0, ldhw, elem(bnt, id0))

    {flmt, id, p, fit} =
      if fit0.ok do
        # run filt_fast once when fit.ok first becomes true
        flmt = if fit.ok, do: flmt, else: filt_fast(sldhw, vldhw, sm, vm, flmt)

        # bins are sorted by volume: only smaller candidates remain interesting
        {Enum.filter(flmt, &(&1 < id0)), id0, p0, fit0}
      else
        flmt = filt_slow(fit.o, vm, flmt)

        if fit.ok do
          {filt_dcol(sm, id0, flmt), id, p, fit}
        else
          if fit0.o > fit.o or (fit0.o == fit.o and id0 < id) do
            {flmt, id0, p0, fit0}
          else
            {flmt, id, p, fit}
          end
        end
      end

    search(ldhw, bnt, sldhw, vldhw, sm, vm, List.delete(flmt, id0), id, p, fit)
  end

  # when all items can fit into some bin, drop bins that cannot hold all items
  # at once, determined via comparing scale, volume and weight
  defp filt_fast(sldhw, {vall, wall}, sm, vm, flmt) do
    Enum.reject(flmt, fn i ->
      {[bl, bd, bh], bw} = elem(sm, i)
      {bv, bwcap} = elem(vm, i)

      Enum.any?(sldhw, fn {[il, id, ih], iw} ->
        il > bl or id > bd or ih > bh or iw > bw
      end) or vall > bv or wall > bwcap
    end)
  end

  # when some items with volume v0 can fit into some bin, drop bins with
  # volume less than v0
  defp filt_slow(v0, vm, flmt) do
    Enum.filter(flmt, fn i ->
      {bv, _} = elem(vm, i)
      bv >= v0
    end)
  end

  # when all items can fit into some bin but bin id0 cannot fit all, drop bins
  # dominated by id0 - they cannot fit all either
  defp filt_dcol(sm, id0, flmt) do
    {[l0, d0, h0], _w0} = elem(sm, id0)

    Enum.reject(flmt, fn i ->
      {[l, d, h], _w} = elem(sm, i)
      l <= l0 + @tol and d <= d0 + @tol and h <= h0 + @tol
    end)
  end
end
