defmodule Robox.Gbp4dTest do
  use ExUnit.Case
  alias Robox.Gbp4d

  @it %{
    oid: [
      1_428_571,
      1_428_571,
      1_428_571,
      1_428_572,
      1_428_572,
      1_428_572,
      1_428_572,
      1_428_572
    ],
    sku: ["A0A0A0", "A0A0A1", "A0A0A1", "A0A0A0", "A0A0A1", "A0A0A1", "A0A0A2", "A0A0A3"],
    l: [2.140000, 7.240000, 7.240000, 2.140000, 7.240000, 7.240000, 6.000000, 4.000000],
    d: [3.580000, 7.240000, 7.240000, 3.580000, 7.240000, 7.240000, 6.000000, 4.000000],
    h: [4.760000, 2.580000, 2.580000, 4.760000, 2.580000, 2.580000, 6.000000, 4.000000],
    w: [243.0000, 110.0000, 110.0000, 243.0000, 110.0000, 110.0000, 235.0000, 258.0000]
  }

  @bn %{
    id: ["K0001", "K0002", "K0003", "K0004", "K0005"],
    l: [06.0000, 10.0000, 09.0000, 10.0000, 22.0000],
    d: [06.0000, 08.0000, 08.0000, 10.0000, 14.0000],
    h: [06.0000, 06.0000, 07.0000, 10.0000, 09.0000],
    w: [600.000, 600.000, 800.000, 800.000, 800.000]
  }

  defp last5_ldhw do
    [@it.l, @it.d, @it.h, @it.w]
    |> Enum.map(&Enum.slice(&1, -5..-1))
    |> Enum.zip()
    |> Enum.map(fn {l, d, h, w} -> {l, d, h, w} end)
  end

  defp bin(i) do
    {Enum.at(@bn.l, i), Enum.at(@bn.d, i), Enum.at(@bn.h, i), Enum.at(@bn.w, i)}
  end

  test "match gbp4d example" do
    ldhw = last5_ldhw()
    m = bin(3)

    p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, m)
    sn = Gbp4d.gbp4d_solver_dpp(p, ldhw, m)

    assert Gbp4d.Ck.gbp4d_checkr(sn)

    assert sn.bn == m
    assert sn.k == [0, 1, 1, 1, 1]
    assert_in_delta sn.o, 550.474816, 1.0e-6
    refute sn.ok

    # heaviest item left out by the weight knapsack, the rest fitted:
    # 2 x (7.24, 7.24, 2.58, 110), (6, 6, 6, 235), (4, 4, 4, 258)
    [it0, it1, it2, it3, it4] = sn.it

    assert {-1.0, -1.0, -1.0, -1.0, 2.14, 3.58, 4.76, 243.0} == it0

    for {x, y, z, _, l, d, h, _} <- [it1, it2, it3, it4] do
      assert x + l <= 10.0 + 1.0e-8
      assert y + d <= 10.0 + 1.0e-8
      assert z + h <= 10.0 + 1.0e-8
    end
  end

  test "match gbp4d example with first bin" do
    ldhw = last5_ldhw()
    m = bin(0)

    p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, m)
    sn = Gbp4d.gbp4d_solver_dpp(p, ldhw, m)

    assert Gbp4d.Ck.gbp4d_checkr(sn)

    # only the 6 x 6 x 6 item fits the 6 x 6 x 6 bin
    assert sn.k == [0, 0, 0, 1, 0]
    assert sn.o == 216.0
    refute sn.ok

    assert Enum.at(sn.it, 3) == {0.0, 0.0, 0.0, 0.0, 6.0, 6.0, 6.0, 235.0}
  end

  test "pack one item in one bin" do
    sn = Gbp4d.pack([[0.1, 0.1, 0.1, 21.8]], [100.0, 100.0, 100.0, 22.0])

    assert sn.k == [1]
    assert sn.ok
    assert_in_delta sn.o, 0.001, 1.0e-9
    assert [{0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 21.8}] == sn.it
  end

  test "pack one item in one bin v2 - item exceeds weight limit" do
    sn = Gbp4d.pack([[0.1, 0.1, 0.1, 22.1]], [100.0, 100.0, 100.0, 22.0])

    assert sn.k == [0]
    refute sn.ok
    assert sn.o == 0.0
    assert [{-1.0, -1.0, -1.0, -1.0, 0.1, 0.1, 0.1, 22.1}] == sn.it
  end

  test "pack one item in one bin v3 - flat item" do
    sn = Gbp4d.pack([[1, 100, 100, 1]], [100.0, 100.0, 100.0, 22.0])

    assert sn.k == [1]
    assert sn.ok
    assert sn.o == 10_000.0
    assert [{0.0, 0.0, 0.0, 0.0, 1.0, 100.0, 100.0, 1.0}] == sn.it
  end

  test "pack multiple item into one bin" do
    sn =
      Gbp4d.pack(
        [
          [40, 30, 30, 8],
          [20, 10, 10, 1],
          [5, 5, 5, 1]
        ],
        [100.0, 100.0, 100.0, 22.0]
      )

    assert Gbp4d.Ck.gbp4d_checkr(sn)
    assert sn.k == [1, 1, 1]
    assert sn.ok
    assert sn.o == 38_125.0

    # exact placement is a score tie broken arbitrarily (the C++ reference
    # stacks along y, this port along z) - assert cumulative weights instead
    assert Enum.map(sn.it, &elem(&1, 3)) == [0.0, 8.0, 9.0]
  end

  test "pack multiple item into one bin v2" do
    sn =
      Gbp4d.pack(
        [
          [60, 40, 40, 10],
          [60, 40, 40, 10],
          [10, 7, 9, 0.44],
          [10, 7, 9, 0.44],
          [10, 7, 9, 0.44],
          [10, 7, 9, 0.44]
        ],
        [100.0, 100.0, 100.0, 22.0]
      )

    assert Gbp4d.Ck.gbp4d_checkr(sn)
    assert sn.k == [1, 1, 1, 1, 1, 1]
    assert sn.ok
    assert_in_delta sn.o, 194_520.0, 1.0e-6
  end

  test "items that cannot all fit keep best subset" do
    # bin volume only fits two of the three 4x4x4 items
    sn = Gbp4d.pack([[4, 4, 4, 1], [4, 4, 4, 1], [4, 4, 4, 1]], [8.0, 4.0, 4.0, 100.0])

    assert Gbp4d.Ck.gbp4d_checkr(sn)
    assert Enum.sum(sn.k) == 2
    assert sn.o == 128.0
    refute sn.ok
  end

  test "zero-dimension item does not crash" do
    # production data occasionally has items with missing dimensions;
    # the old Matrex port crashed on these (NaN in the profit ranking)
    sn = Gbp4d.pack([[0, 0, 0, 1], [4, 4, 4, 1]], [8.0, 4.0, 4.0, 100.0])

    assert Gbp4d.Ck.gbp4d_checkr(sn)
    # matches the C++ reference: a zero-volume item never raises the
    # objective, so it is reported as not fitted
    assert sn.k == [0, 1]
    assert sn.o == 64.0
  end

  test "empty item list" do
    sn = Gbp4d.pack([], [10.0, 10.0, 10.0, 100.0])

    assert sn.k == []
    assert sn.o == 0.0
    assert sn.ok
  end
end
