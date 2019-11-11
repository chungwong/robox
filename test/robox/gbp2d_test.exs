defmodule Robox.Gbp2dTest do
  use ExUnit.Case
  alias Robox.{Gbp2d}

  test "match gbp2d example" do
    ldhw =
      Matrex.new([
        [2.14,  7.24,  7.24, 6,  4],
        [3.58,  7.24,  7.24, 6,  4],
      ])

    m =
      Matrex.new([
        [10],
        [10],
      ])

    p = Matrex.new([
      [0],
      [3],
      [4],
      [2],
      [1]
    ])

    sn = Gbp2d.gbp2d_solver_dpp(p, ldhw, m)

    assert %Gbp2d{
             bn: m,
             it:
               Matrex.new([
                 [7.24, -1.0,  0.0, -1.0,  -1],
                 [0.0,  -1.0,  0.0, -1.0,  -1],
                 [2.14, 7.24, 7.24,  6.0, 4.0],
                 [3.58, 7.24, 7.24,  6.0, 4.0],
               ]),
             k:
               Matrex.new([
                 [1.0],
                 [0.0],
                 [1.0],
                 [0.0],
                 [0.0]
               ]),
             o: 60.07879638671875,
             ok: false,
             p: p
           } == sn
  end

  test "2d: 1 item 1 bin" do
    ldhw =
      Matrex.new([
        [0.1],
        [0.1],
      ])

    m =
      Matrex.new([
        [10],
        [10],
      ])

    p = Matrex.new([
      [0],
    ])

    sn = Gbp2d.gbp2d_solver_dpp(p, ldhw, m)

    assert %Gbp2d{
             bn: m,
             it:
               Matrex.new([
                 [0],
                 [0.0],
                 [0.1],
                 [0.1],
               ]),
             k:
               Matrex.new([
                 [1.0],
               ]),
             o: 0.010000000707805157,
             ok: true,
             p: p
           } == sn
  end

  test "2d: 1 item 1 bin (too large)" do
    ldhw =
      Matrex.new([
        [0.1],
        [10.1],
      ])

    m =
      Matrex.new([
        [10],
        [10],
      ])

    p = Matrex.new([
      [0],
    ])

    sn = Gbp2d.gbp2d_solver_dpp(p, ldhw, m)

    assert %Gbp2d{
             bn: m,
             it:
               Matrex.new([
                 [-1],
                 [-1],
                 [0.1],
                 [10.1],
               ]),
             k:
               Matrex.new([
                 [0.0],
               ]),
             o: 0.0,
             ok: false,
             p: p
           } == sn
  end


  def gen_products() do
    max_length = 100
    max_depth = 100

    Stream.flat_map(1..trunc(max_length), fn l ->
      Stream.map(1..trunc(max_depth), fn d ->
        [
          l,
          d,
        ]
      end)
    end)
  end

  test "2d: crazy test (1 x 1 to 100 x 100)" do
    m =
      Matrex.new([
        [100],
        [100],
      ])

    p = Matrex.new([
      [0],
    ])

    Enum.each(gen_products(), fn [l, d] ->
      ldhw =
        Matrex.new([
          [l],
          [d],
        ])
      IO.inspect ldhw, label: "crazy test ldhw ="
      sn = Gbp2d.gbp2d_solver_dpp(p, ldhw, m)
      IO.inspect sn
    end)
  end
end

