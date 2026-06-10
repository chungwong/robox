defmodule Robox.Gbp4d.It do
  @moduledoc """
  Candidate placement (ktlist) creation for gbp4d, mirroring gbp/src/gbp4d_it.cpp.
  """

  alias Robox.Gbp4d.{Ktlist4d, Xp}

  @doc """
  Select all feasible placements of the item `ktinit` w.r.t orientation and
  extreme point, score each by the entropy of the residual spaces it leaves,
  and keep the best `nlmt` (all when `nlmt` is 0 or exceeds the count).

  `its` are the items already in the bin.
  """
  @spec gbp4d_it_create_ktlist(tuple, [tuple], [tuple], tuple, non_neg_integer) :: %Ktlist4d{}
  def gbp4d_it_create_ktlist(bn, its, xp, ktinit, nlmt) do
    ktldht = gbp4d_it_create_ktldht(ktinit)

    # fast filt w.r.t extreme point residual space
    kts =
      for {l, d, h, w} <- ktldht,
          {xpx, xpy, xpz, xpw, xpl, xpd, xph, xpwt} <- xp,
          l <= xpl and d <= xpd and h <= xph and w <= xpwt do
        {xpx, xpy, xpz, xpw, l, d, h, w}
      end

    # slow filt w.r.t no conflict with any existing it
    kts =
      Enum.filter(kts, fn {kx, ky, kz, _, kl, kd, kh, _} ->
        Enum.all?(its, fn {ix, iy, iz, _, il, id, ih, _} ->
          kx + kl <= ix or ix + il <= kx or
            ky + kd <= iy or iy + id <= ky or
            kz + kh <= iz or iz + ih <= kz
        end)
      end)

    scored =
      Enum.map(kts, fn kt ->
        xp_with_kt = Xp.gbp4d_xp_update_xp(bn, its, kt, xp)
        {kt, xp_with_kt, gbp4d_it_scorer_ktlist(xp_with_kt)}
      end)

    purified = gbp4d_it_purify_ktlist(scored, nlmt)

    %Ktlist4d{
      n: length(purified),
      kt: Enum.map(purified, &elem(&1, 0)),
      xp: Enum.map(purified, &elem(&1, 1)),
      s: Enum.map(purified, &elem(&1, 2))
    }
  end

  @doc """
  Sort scored placements by score ascending and keep the best `nlmt`
  (all when `nlmt` is 0 or >= the number of placements).
  """
  @spec gbp4d_it_purify_ktlist([{tuple, [tuple], float}], non_neg_integer) ::
          [{tuple, [tuple], float}]
  def gbp4d_it_purify_ktlist(scored, nlmt) do
    sorted = Enum.sort_by(scored, &elem(&1, 2))

    if nlmt == 0 or nlmt >= length(sorted) do
      sorted
    else
      Enum.take(sorted, nlmt)
    end
  end

  @doc """
  Score a placement via the entropy of the residual space of the extreme
  points it leaves - the smaller the better: prefer configurations where few
  large available spaces dominate.
  """
  @spec gbp4d_it_scorer_ktlist([tuple]) :: float
  def gbp4d_it_scorer_ktlist(xp_with_kt) do
    rs = Enum.map(xp_with_kt, fn {_, _, _, _, l, d, h, _} -> l * d * h end)
    total = Enum.sum(rs)

    Enum.reduce(rs, 0.0, fn r, acc ->
      r = r / total
      acc - r * :math.log(r)
    end)
  end

  @doc """
  Create all possible non-duplicated orientations of `ktinit` as a list of
  `{l, d, h, w}` tuples.
  """
  @spec gbp4d_it_create_ktldht(tuple) :: [tuple]
  def gbp4d_it_create_ktldht({_, _, _, _, l, d, h, w}) do
    all = [
      {l, d, h, w},
      {l, h, d, w},
      {d, l, h, w},
      {d, h, l, w},
      {h, l, d, w},
      {h, d, l, w}
    ]

    indices =
      cond do
        l == d and d == h -> [0]
        l == d -> [0, 1, 4]
        d == h -> [0, 2, 3]
        h == l -> [0, 1, 2]
        true -> [0, 1, 2, 3, 4, 5]
      end

    Enum.map(indices, &Enum.at(all, &1))
  end
end
