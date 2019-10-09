defmodule Robox.Gbp1dTest do
  use ExUnit.Case
  alias Robox.{Gbp1d}

  test "match gbp1d example" do
    p = Matrex.new([[36.46731, 135.23741, 135.23741, 216.00000, 64.00000]])
    w = Matrex.new([[243, 110, 100, 235, 258]])
    c = 714.28

    assert %Gbp1d{
             c: trunc(c),
             k:
               Matrex.new([
                 [0.0],
                 [1.0],
                 [1.0],
                 [1.0],
                 [1.0]
               ]),
             o: 550.474853515625,
             ok: false,
             w: w,
             p: p
           } == Gbp1d.gbp1d_solver_dpp(p, w, c)
  end
end
