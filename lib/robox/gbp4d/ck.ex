defmodule Robox.Gbp4d.Ck do
  @moduledoc """
  Solution validation for gbp4d / gbp4q, mirroring gbp/src/gbp4d_ck.cpp.
  """

  require Logger

  @doc """
  Check a `Robox.Gbp4d` solution is valid: no conflict between item and bin,
  no conflict between each pair of items, and no conflict on weight limit.
  """
  @spec gbp4d_checkr(%Robox.Gbp4d{}) :: boolean
  def gbp4d_checkr(%{it: it, bn: {bl, bd, bh, bw}, k: k}) do
    fitted =
      Enum.zip(it, k)
      |> Enum.filter(fn {_, ki} -> ki == 1 end)
      |> Enum.map(&elem(&1, 0))

    ok_bn =
      Enum.all?(fitted, fn {x, y, z, _, l, d, h, _} ->
        within = x + l <= bl and y + d <= bd and z + h <= bh

        unless within do
          Logger.warning("gbp4d_checkr: it conflict bn: #{inspect({x, y, z, l, d, h})}")
        end

        within
      end)

    ok_it =
      ok_bn and
        no_pairwise_conflict?(fitted)

    weight = Enum.reduce(fitted, 0.0, fn {_, _, _, _, _, _, _, w}, acc -> acc + w end)

    ok_w = weight <= bw

    unless ok_w do
      Logger.warning("gbp4d_checkr: it conflict bn: conflict on weight constraint.")
    end

    ok_it and ok_w
  end

  @doc """
  Check a `Robox.Gbp4q` solution is valid w.r.t the single selected bin.
  """
  @spec gbp4q_checkr(%Robox.Gbp4q{}) :: boolean
  def gbp4q_checkr(%Robox.Gbp4q{} = sn) do
    case Enum.with_index(sn.f) |> Enum.filter(fn {fi, _} -> fi == 1 end) do
      [{_, id}] ->
        {bl, bd, bh, bw} = Enum.at(sn.bn, id)

        gbp4d_checkr(%{it: sn.it, bn: {bl, bd, bh, bw}, k: sn.k})

      _ ->
        Logger.warning("gbp4q_checkr: f should have a unique index label 1.")
        false
    end
  end

  defp no_pairwise_conflict?([]), do: true

  defp no_pairwise_conflict?([{ix, iy, iz, _, il, id, ih, _} = _it | rest]) do
    ok =
      Enum.all?(rest, fn {jx, jy, jz, _, jl, jd, jh, _} ->
        sep =
          ix + il <= jx or jx + jl <= ix or
            iy + id <= jy or jy + jd <= iy or
            iz + ih <= jz or jz + jh <= iz

        unless sep do
          Logger.warning("gbp4d_checkr: it conflict it.")
        end

        sep
      end)

    ok and no_pairwise_conflict?(rest)
  end
end
