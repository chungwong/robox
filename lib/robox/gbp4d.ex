defmodule Robox.Gbp4d do
  @moduledoc """
  Generalized bin packing problem in 4 dimensions, a.k.a bin packing problem
  with weight limit, mirroring gbp/src/gbp4d.cpp.

  Solves

      maximize   sum_{j=1}^{n} p_j k_j

      subject to sum_{j=1}^{n} w_j k_j <= mw and
                 fit (l_j, d_j, h_j) at coordinate (x_j, y_j, z_j)
                 such that no overlap in ml x md x mh cuboid

  via extreme point heuristic and best information score fit strategy.

  Result struct:

    * `p`  - profit vector used as fit sequence ranking, list of floats
    * `it` - one entry per item, 8-tuple `{x, y, z, w, l, d, h, wt}` where
      `x, y, z` is the fitted position (`-1.0` when not fitted) and `w` the
      cumulative weight in the bin when the item was placed
    * `bn` - bin as `{l, d, h, w}` tuple
    * `k`  - selection indicator per item, 0 | 1
    * `o`  - objective: total volume of fitted items
    * `ok` - did all items fit?
  """

  alias Robox.Gbp1d
  alias Robox.Gbp4d.It

  defstruct [:p, :it, :bn, :k, :o, :ok]

  @type t :: %__MODULE__{}

  # gbp4d global setup on parameter
  # when <  2 items search over all in ktlist w.r.t orientation + extreme point
  @gbp4d_nlvl_mkt0 2
  # when <  3 items search best   5 in ktlist w.r.t orientation + extreme point
  @gbp4d_nlvl_mkt1 3
  # when <  5 items search best   3 in ktlist w.r.t orientation + extreme point
  @gbp4d_nlvl_mkt2 5
  # when <  8 items search best   2 in ktlist w.r.t orientation + extreme point
  @gbp4d_nlvl_mkt3 8
  # when >= 8 items search best   1 in ktlist w.r.t orientation + extreme point

  @doc """
  Pack items `ldhw` (list of `[l, d, h, w]` per item) into the single bin
  `bn` (`[l, d, h, w]` or `[[l, d, h, w]]`).
  """
  @spec pack([[number]], [[number]] | [number]) :: t
  def pack(ldhw, [bn]) when is_list(bn), do: pack(ldhw, bn)

  def pack(ldhw, [l, d, h, w]) do
    ldhw = Enum.map(ldhw, fn [il, id, ih, iw] -> {il / 1, id / 1, ih / 1, iw / 1} end)
    m = {l / 1, d / 1, h / 1, w / 1}

    p = gbp4d_solver_dpp_prep_create_p(ldhw, m)
    gbp4d_solver_dpp(p, ldhw, m)
  end

  @doc """
  Solve gbp4d via extreme point heuristic and best information score fit
  strategy. `p` is the profit list, `ldhw` a list of `{l, d, h, w}` tuples,
  `m` the bin as `{l, d, h, w}`.
  """
  @spec gbp4d_solver_dpp([float], [tuple], tuple) :: t
  def gbp4d_solver_dpp(p, ldhw, {ml, md, mh, mw}) do
    n = length(p)
    bn = {ml / 1, md / 1, mh / 1, mw / 1}

    # it: id => {x, y, z, w, l, d, h, wt}, position initialised at origin
    it =
      ldhw
      |> Enum.with_index()
      |> Map.new(fn {{l, d, h, w}, i} -> {i, {0.0, 0.0, 0.0, 0.0, l / 1, d / 1, h / 1, w / 1}} end)

    if n == 0 do
      %__MODULE__{p: p, it: [], bn: bn, k: [], o: 0.0, ok: true}
    else
      # q: it fit sequence - descending profit
      q =
        0..(n - 1)
        |> Enum.sort_by(fn i -> -Enum.at(p, i) end)
        |> List.to_tuple()

      xp = [{0.0, 0.0, 0.0, 0.0, ml / 1, md / 1, mh / 1, mw / 1}]

      # nastr: number of it certain to be left over by volume or weight limit
      nastr = gbp4d_solver_dpp_main_create_nastr(p, ldhw, bn)

      # g: status per id: 0 open, 1 fitted, 2 determined no fit
      g = Map.new(0..(n - 1), &{&1, 0})

      # v: volume per id - l, d, h in objective and constraint, w in constraint only
      v = Map.new(it, fn {i, {_, _, _, _, l, d, h, _}} -> {i, l * d * h} end)

      # vastr: maximum volume achievable via fit
      vastr = min(Enum.sum(Map.values(v)), ml * md * mh)

      ctx = %{bn: bn, q: q, n: n, nastr: nastr, v: v, vastr: vastr}

      # global best over the recursive search
      best = %{itastr: it, gastr: g, uastr: 0.0}

      {_ok, best} = gbp4d_solver_dpp_main(ctx, it, xp, n, g, 0.0, best)

      k = Enum.map(0..(n - 1), fn i -> if best.gastr[i] == 1, do: 1, else: 0 end)

      o =
        Enum.reduce(0..(n - 1), 0.0, fn i, acc ->
          if best.gastr[i] == 1, do: acc + v[i], else: acc
        end)

      ok = Enum.all?(k, &(&1 == 1))

      # flag no fit with (x, y, z, w) = (-1, -1, -1, -1) instead of (0, 0, 0, 0)
      it_out =
        Enum.map(0..(n - 1), fn i ->
          if best.gastr[i] == 1 do
            best.itastr[i]
          else
            {_, _, _, _, l, d, h, w} = it[i]
            {-1.0, -1.0, -1.0, -1.0, l, d, h, w}
          end
        end)

      %__MODULE__{p: p, it: it_out, bn: bn, k: k, o: o, ok: ok}
    end
  end

  @doc """
  Recursive main solver. Returns `{ok, best}` where `ok` means all items fit
  into the bin or a subset of items makes full use of the bin volume, and
  `best` tracks the global best `itastr` / `gastr` / `uastr`.
  """
  def gbp4d_solver_dpp_main(ctx, it, xp, nlvl, g, u, best) do
    # id: index of it to fit at this level
    id = elem(ctx.q, ctx.n - nlvl)

    # glmt: items neither fitted nor determined no fit at entry - for vmiss
    glmt = for {i, 0} <- g, do: i

    it0 = for {i, 1} <- g, do: it[i]
    kt0 = it[id]
    nlmt = gbp4d_solver_dpp_main_create_nlmt(nlvl, ctx.nastr)

    ktlist = It.gbp4d_it_create_ktlist(ctx.bn, it0, xp, kt0, nlmt)

    cond do
      nlvl == 1 ->
        if ktlist.n > 0 do
          it = Map.put(it, id, hd(ktlist.kt))
          g = Map.put(g, id, 1)
          u = u + ctx.v[id]

          best = if u > best.uastr, do: %{itastr: it, gastr: g, uastr: u}, else: best

          {best.uastr == ctx.vastr or all_fitted?(best.gastr), best}
        else
          {false, best}
        end

      ktlist.n > 0 ->
        g1 = Map.put(g, id, 1)
        u1 = u + ctx.v[id]

        {ok, best} =
          Enum.zip(ktlist.kt, ktlist.xp)
          |> Enum.reduce_while({false, best}, fn {kt, xp1}, {_ok, best} ->
            it1 = Map.put(it, id, kt)

            best = if u1 > best.uastr, do: %{itastr: it1, gastr: g1, uastr: u1}, else: best

            if best.uastr == ctx.vastr or all_fitted?(best.gastr) do
              {:halt, {true, best}}
            else
              case gbp4d_solver_dpp_main(ctx, it1, xp1, nlvl - 1, g1, u1, best) do
                {true, best} -> {:halt, {true, best}}
                {false, best} -> {:cont, {false, best}}
              end
            end
          end)

        # what if skip this it? can later it combined be better?
        if not ok and nlvl < @gbp4d_nlvl_mkt2 do
          vmiss =
            glmt
            |> Enum.filter(fn i -> best.gastr[i] != 1 end)
            |> Enum.reduce(0.0, fn i, acc -> acc + ctx.v[i] end)

          if vmiss > ctx.v[id] do
            g2 = Map.put(g, id, 2)
            gbp4d_solver_dpp_main(ctx, it, xp, nlvl - 1, g2, u, best)
          else
            {false, best}
          end
        else
          {ok, best}
        end

      true ->
        g2 = Map.put(g, id, 2)
        gbp4d_solver_dpp_main(ctx, it, xp, nlvl - 1, g2, u, best)
    end
  end

  defp all_fitted?(gastr) do
    Enum.all?(gastr, fn {_, s} -> s == 1 end)
  end

  @doc """
  Number of items certain to be left over due to the bin volume or weight
  limit, walking the fit sequence accumulating volume and weight.
  """
  @spec gbp4d_solver_dpp_main_create_nastr([float], [tuple], tuple) :: non_neg_integer
  def gbp4d_solver_dpp_main_create_nastr(p, ldhw, {ml, md, mh, mw}) do
    n = length(p)

    q =
      0..(n - 1)
      |> Enum.sort_by(fn i -> -Enum.at(p, i) end)

    v = Enum.map(ldhw, fn {l, d, h, _} -> l * d * h end) |> List.to_tuple()
    w = Enum.map(ldhw, fn {_, _, _, wt} -> wt / 1 end) |> List.to_tuple()

    mv = ml * md * mh

    {nastr, _, _, _} =
      Enum.reduce_while(Enum.with_index(q), {0, 0.0, 0.0, n}, fn {qi, i}, {nastr, v0, w0, n} ->
        v0 = v0 + elem(v, qi)
        w0 = w0 + elem(w, qi)

        if v0 >= mv or w0 >= mw do
          {:halt, {n - 1 - i, v0, w0, n}}
        else
          {:cont, {nastr, v0, w0, n}}
        end
      end)

    nastr
  end

  @doc """
  Limit on the number of candidate placements searched per level.
  """
  @spec gbp4d_solver_dpp_main_create_nlmt(non_neg_integer, non_neg_integer) :: non_neg_integer
  def gbp4d_solver_dpp_main_create_nlmt(nlvl, nastr) do
    cond do
      nlvl <= nastr -> 1
      nlvl - nastr < @gbp4d_nlvl_mkt0 -> 0
      nlvl - nastr < @gbp4d_nlvl_mkt1 -> 5
      nlvl - nastr < @gbp4d_nlvl_mkt2 -> 3
      nlvl - nastr < @gbp4d_nlvl_mkt3 -> 2
      true -> 1
    end
  end

  @doc """
  Create the profit vector `p` via clustering on weight (knapsack fit),
  max(l, d, h) bucket and area, so that items with higher `p` get fit into
  the bin earlier in the queue.

  `ldhw` is a list of `{l, d, h, w}` tuples, `m` the bin `{l, d, h, w}`.
  """
  @spec gbp4d_solver_dpp_prep_create_p([tuple], tuple) :: [float]
  def gbp4d_solver_dpp_prep_create_p(ldhw, m)

  def gbp4d_solver_dpp_prep_create_p([], _m), do: []

  def gbp4d_solver_dpp_prep_create_p(ldhw, {ml, md, mh, mw}) do
    # q0: cluster on weight w.r.t gbp1d knapsack at 0.25 granularity
    v = Enum.map(ldhw, fn {l, d, h, _} -> l * d * h end)
    w = Enum.map(ldhw, fn {_, _, _, wt} -> ceil(wt / 0.25) end)
    c = trunc(Float.floor(mw / 0.25))

    wfit = Gbp1d.gbp1d_solver_dpp(v, w, c)

    q0 = Enum.map(wfit.k, fn ki -> if ki == 1, do: 1.0, else: 0.0 end)

    # q1: cluster on max(l, d, h) w.r.t bn max(l, d, h) into per 0.25 inch
    mldh = Enum.max([ml, md, mh])
    ildh = Enum.map(ldhw, fn {l, d, h, _} -> Enum.max([l, d, h]) end)

    nldh = linspace(0.0, mldh, ceil(mldh * 4.0))

    q1 =
      Enum.map(ildh, fn iv ->
        # sum of 0-based indices of nldh values <= iv
        k = Enum.count(nldh, &(&1 <= iv))
        k * (k - 1) / 2
      end)

    # q2: area - volume over max dimension (0-dimension items rank lowest)
    q2 = Enum.zip(v, ildh) |> Enum.map(fn {vi, ii} -> if ii == 0.0, do: 0.0, else: vi / ii end)

    # rank items by (q0, q1, q2) ascending - higher rank fits earlier
    order =
      Enum.zip([q0, q1, q2])
      |> Enum.with_index()
      |> Enum.sort_by(fn {key, _i} -> key end)
      |> Enum.map(fn {_key, i} -> i end)

    order
    |> Enum.with_index()
    |> Enum.sort_by(fn {i, _rank} -> i end)
    |> Enum.map(fn {_i, rank} -> rank / 1 end)
  end

  # arma::linspace - num evenly spaced values from start to stop inclusive
  defp linspace(_start, stop, num) when num <= 1, do: [stop]

  defp linspace(start, stop, num) do
    step = (stop - start) / (num - 1)
    Enum.map(0..(num - 1), fn i -> start + step * i end)
  end
end
