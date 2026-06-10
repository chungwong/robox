defmodule Robox.Gbp1dTest do
  use ExUnit.Case
  alias Robox.Gbp1d

  test "match gbp1d example" do
    p = [36.46731, 135.23741, 135.23741, 216.0, 64.0]
    w = [243, 110, 100, 235, 258]
    c = 714.28

    sn = Gbp1d.gbp1d_solver_dpp(p, w, c)

    assert sn.c == 714
    assert sn.k == [0, 1, 1, 1, 1]
    assert_in_delta sn.o, 550.47482, 1.0e-8
    refute sn.ok
  end

  test "all items fit" do
    sn = Gbp1d.gbp1d_solver_dpp([1.0, 2.0], [3, 4], 10)

    assert sn.k == [1, 1]
    assert sn.o == 3.0
    assert sn.ok
  end

  test "no items" do
    sn = Gbp1d.gbp1d_solver_dpp([], [], 10)

    assert sn.k == []
    assert sn.o == 0.0
    assert sn.ok
  end

  test "item heavier than capacity" do
    sn = Gbp1d.gbp1d_solver_dpp([5.0], [20], 10)

    assert sn.k == [0]
    assert sn.o == 0.0
    refute sn.ok
  end
end
