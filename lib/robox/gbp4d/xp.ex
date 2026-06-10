defmodule Robox.Gbp4d.Xp do
  @moduledoc """
  Extreme point calculation for gbp4d, mirroring gbp/src/gbp4d_xp.cpp.

  An item `it`, a candidate `kt` and an extreme point `xp` are all 8-element
  tuples `{x, y, z, w, l, d, h, wt}`:

    * `x, y, z` position in the bin, `w` cumulative weight held when placed
    * `l, d, h` scale along x, y, z; `wt` weight (for an extreme point these
      are the residual space along x, y, z and the residual weight)

  A bin `bn` is a 4-element tuple `{l, d, h, w}`.
  """

  @tol 1.0e-8

  @doc """
  Calculate the extreme points of a bin holding the items `its` from scratch.
  """
  @spec gbp4d_xp_create_xp(tuple, [tuple]) :: [tuple]
  def gbp4d_xp_create_xp({bl, bd, bh, bw}, its) do
    # fit items one by one in z, y, x order - mimic fit sequence
    sorted = Enum.sort_by(its, fn {x, y, z, _, _, _, _, _} -> {z, y, x} end)

    xp0 = [{0.0, 0.0, 0.0, 0.0, bl / 1, bd / 1, bh / 1, bw / 1}]

    sorted
    |> Enum.with_index()
    |> Enum.reduce(xp0, fn {kt, i}, xp ->
      gbp4d_xp_update_xp({bl, bd, bh, bw}, Enum.take(sorted, i), kt, xp)
    end)
  end

  @doc """
  Update the extreme point list `xp` after fitting item `kt` into the bin `bn`
  already holding items `its` (`its` must not include `kt`).
  """
  @spec gbp4d_xp_update_xp(tuple, [tuple], tuple, [tuple]) :: [tuple]
  def gbp4d_xp_update_xp({bl, bd, bh, bw}, its, kt, xp) do
    # remove extreme points taken by or fallen into kt and
    # update residual space of surviving extreme points w.r.t kt
    xp = gbp4d_xp_update_xp_ikt(kt, xp)

    {ktx, kty, ktz, ktw, ktl, ktd, kth, ktwt} = kt

    # maxBound: 6 projection positions, from the 3 corner points of kt
    # projected back onto items (or the bin walls at 0)
    mb = Enum.reduce(its, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, &gbp4d_xp_update_maxbnd(&1, kt, &2))
    {mb0, mb1, mb2, mb3, mb4, mb5} = mb

    wp = ktw + ktwt

    # 6 candidate new extreme points
    positions = [
      {ktx + ktl, mb0, ktz, wp},
      {ktx + ktl, kty, mb1, wp},
      {ktx, kty + ktd, mb2, wp},
      {mb3, kty + ktd, ktz, wp},
      {mb4, kty, ktz + kth, wp},
      {ktx, mb5, ktz + kth, wp}
    ]

    # minBound: residual space upper corner for each candidate, shrunk by
    # items standing between the candidate position and the bin walls
    min_bounds =
      Enum.map(positions, fn pos ->
        Enum.reduce(its, {bl / 1, bd / 1, bh / 1, bw / 1}, fn it, minb ->
          gbp4d_xp_update_minbnd(it, pos, minb)
        end)
      end)

    candidates =
      Enum.zip(positions, min_bounds)
      |> Enum.map(fn {{px, py, pz, pw}, {mx, my, mz, mw}} ->
        {px, py, pz, pw, mx - px, my - py, mz - pz, mw - pw}
      end)

    (xp ++ candidates)
    |> unique_cols()
    |> gbp4d_xp_purify_xp()
    |> Enum.sort_by(fn {x, y, z, _, _, _, _, _} -> {z, y, x} end)
  end

  @doc """
  Remove extreme points with zero residual space and extreme points dominated
  by another extreme point in the list.
  """
  @spec gbp4d_xp_purify_xp([tuple]) :: [tuple]
  def gbp4d_xp_purify_xp(xp) do
    xp =
      Enum.reject(xp, fn {_, _, _, _, l, d, h, w} ->
        l == 0.0 or d == 0.0 or h == 0.0 or w == 0.0
      end)

    Enum.reject(xp, fn {xj, yj, zj, _, lj, dj, hj, _} = j ->
      Enum.any?(xp, fn {xi, yi, zi, _, li, di, hi, _} = i ->
        i != j and xi <= xj and yi <= yj and zi <= zj and
          li >= lj and di >= dj and hi >= hj
      end)
    end)
  end

  @doc """
  Remove extreme points taken by or fallen into `kt`, then update the residual
  space of the remaining extreme points w.r.t `kt`.
  """
  @spec gbp4d_xp_update_xp_ikt(tuple, [tuple]) :: [tuple]
  def gbp4d_xp_update_xp_ikt({ktx, kty, ktz, _ktw, ktl, ktd, kth, ktwt}, xp) do
    xp
    |> Enum.reject(fn {x, y, z, _, _, _, _, _} ->
      ktx <= x and x < ktx + ktl and
        kty <= y and y < kty + ktd and
        ktz <= z and z < ktz + kth
    end)
    |> Enum.map(fn {x, y, z, w, l, d, h, wt} ->
      l =
        if x <= ktx and y >= kty and y < kty + ktd and z >= ktz and z < ktz + kth do
          min(l, ktx - x)
        else
          l
        end

      d =
        if y <= kty and z >= ktz and z < ktz + kth and x >= ktx and x < ktx + ktl do
          min(d, kty - y)
        else
          d
        end

      h =
        if z <= ktz and x >= ktx and x < ktx + ktl and y >= kty and y < kty + ktd do
          min(h, ktz - z)
        else
          h
        end

      # weight on separate single dimension: weight holding grows, available shrinks
      {x, y, z, w + ktwt, l, d, h, wt - ktwt}
    end)
  end

  @doc """
  Shrink the residual space bound `{mx, my, mz, mw}` of a candidate extreme
  point at `pos` w.r.t a single item `it` standing in the way.
  """
  @spec gbp4d_xp_update_minbnd(tuple, tuple, tuple) :: tuple
  def gbp4d_xp_update_minbnd(it, {px, py, pz, pw}, {mx, my, mz, mw}) do
    {itx, ity, itz, _, _, _, _, _} = it

    # a virtual kt as a single point with no scale at the candidate position
    akt = {px, py, pz, pw, 0.0, 0.0, 0.0, 0.0}

    {ik0, ik1, ik2, ik3, ik4, ik5} = gbp4d_xp_it_qjt_kt(it, akt)

    mx = if ik3 and ik4, do: min(itx, mx), else: mx
    my = if ik5 and ik0, do: min(ity, my), else: my
    mz = if ik1 and ik2, do: min(itz, mz), else: mz

    {mx, my, mz, mw}
  end

  @doc """
  Raise the projection bounds `maxBound` of the new extreme points spawned by
  `kt` w.r.t a single item `it` that can take the projection.
  """
  @spec gbp4d_xp_update_maxbnd(tuple, tuple, tuple) :: tuple
  def gbp4d_xp_update_maxbnd(it, kt, {mb0, mb1, mb2, mb3, mb4, mb5}) do
    {itx, ity, itz, _, itl, itd, ith, _} = it

    {ik0, ik1, ik2, ik3, ik4, ik5} = gbp4d_xp_it_pjt_kt(it, kt)

    mb0 = if ik0 and ity + itd > mb0, do: ity + itd, else: mb0
    mb1 = if ik1 and itz + ith > mb1, do: itz + ith, else: mb1
    mb2 = if ik2 and itz + ith > mb2, do: itz + ith, else: mb2
    mb3 = if ik3 and itx + itl > mb3, do: itx + itl, else: mb3
    mb4 = if ik4 and itx + itl > mb4, do: itx + itl, else: mb4
    mb5 = if ik5 and ity + itd > mb5, do: ity + itd, else: mb5

    {mb0, mb1, mb2, mb3, mb4, mb5}
  end

  @doc """
  Can item `it` take the projection of the extreme points created by `kt`, on
  the reverse direction (directions xY, xZ, yZ, yX, zX, zY).
  """
  @spec gbp4d_xp_it_qjt_kt(tuple, tuple) :: tuple
  def gbp4d_xp_it_qjt_kt(
        {itx, ity, itz, _, itl, itd, ith, _},
        {ktx, kty, ktz, _, ktl, ktd, kth, _}
      ) do
    {
      kty + ktd <= ity and itx <= ktx + ktl and ktx + ktl < itx + itl and
        itz <= ktz and ktz < itz + ith,
      ktz + kth <= itz and itx <= ktx + ktl and ktx + ktl < itx + itl and
        ity <= kty and kty < ity + itd,
      ktz + kth <= itz and ity <= kty + ktd and kty + ktd < ity + itd and
        itx <= ktx and ktx < itx + itl,
      ktx + ktl <= itx and ity <= kty + ktd and kty + ktd < ity + itd and
        itz <= ktz and ktz < itz + ith,
      ktx + ktl <= itx and itz <= ktz + kth and ktz + kth < itz + ith and
        ity <= kty and kty < ity + itd,
      kty + ktd <= ity and itz <= ktz + kth and ktz + kth < itz + ith and
        itx <= ktx and ktx < itx + itl
    }
  end

  @doc """
  Can item `it` take the projection of the extreme points created by `kt`, on
  the one direction (directions xY, xZ, yZ, yX, zX, zY).
  """
  @spec gbp4d_xp_it_pjt_kt(tuple, tuple) :: tuple
  def gbp4d_xp_it_pjt_kt(
        {itx, ity, itz, _, itl, itd, ith, _},
        {ktx, kty, ktz, _, ktl, ktd, kth, _}
      ) do
    {
      ity + itd <= kty and itx <= ktx + ktl and ktx + ktl < itx + itl and
        itz <= ktz and ktz < itz + ith,
      itz + ith <= ktz and itx <= ktx + ktl and ktx + ktl < itx + itl and
        ity <= kty and kty < ity + itd,
      itz + ith <= ktz and ity <= kty + ktd and kty + ktd < ity + itd and
        itx <= ktx and ktx < itx + itl,
      itx + itl <= ktx and ity <= kty + ktd and kty + ktd < ity + itd and
        itz <= ktz and ktz < itz + ith,
      itx + itl <= ktx and itz <= ktz + kth and ktz + kth < itz + ith and
        ity <= kty and kty < ity + itd,
      ity + itd <= kty and itz <= ktz + kth and ktz + kth < itz + ith and
        itx <= ktx and ktx < itx + itl
    }
  end

  @doc """
  Drop later entries approximately equal (absdiff within tolerance) to an
  earlier entry, keeping first occurrences - mirrors `unique_cols`.
  """
  @spec unique_cols([tuple]) :: [tuple]
  def unique_cols(cols) do
    Enum.reduce(cols, [], fn col, kept ->
      if Enum.any?(kept, &approx_equal?(&1, col)) do
        kept
      else
        [col | kept]
      end
    end)
    |> Enum.reverse()
  end

  defp approx_equal?(a, b) when tuple_size(a) == tuple_size(b) do
    Enum.all?(0..(tuple_size(a) - 1), fn i -> abs(elem(a, i) - elem(b, i)) <= @tol end)
  end

  defp approx_equal?(_, _), do: false
end
