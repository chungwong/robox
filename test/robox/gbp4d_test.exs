defmodule Robox.Gbp4dTest do
  use ExUnit.Case
  alias Robox.{Gbp4d}

  require Logger

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

  test "match gbp4d example" do
    ldhw =
      Matrex.new([
        @it.l |> Enum.slice(-5..-1),
        @it.d |> Enum.slice(-5..-1),
        @it.h |> Enum.slice(-5..-1),
        @it.w |> Enum.slice(-5..-1)
      ])

    m =
      Matrex.new([
        @bn.l,
        @bn.d,
        @bn.h,
        @bn.w
      ])

    m4l = Matrex.column(m, 4)

    p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, m4l)

    sn = Gbp4d.gbp4d_solver_dpp(p, ldhw, m4l)
    Gbp4d.Ck.gbp4d_checkr(sn)

    assert %Gbp4d{
             bn: m4l,
             it:
               Matrex.new([
                 [-1.0, 0.0, 0.0, 0.0, 6.0],
                 [-1.0, 0.0, 0.0, 2.58, 2.58],
                 [-1.0, 2.58, 0.0, 2.58, 2.58],
                 [-1.0, 110.0, 0.0, 220.0, 455.0],
                 [2.14, 7.24, 7.24, 6.0, 4.0],
                 [4.76, 2.58, 7.24, 6.0, 4.0],
                 [3.58, 7.24, 2.58, 6.0, 4.0],
                 [243.0, 110.0, 110.0, 235.0, 258.0]
               ]),
             k:
               Matrex.new([
                 [0.0],
                 [1.0],
                 [1.0],
                 [1.0],
                 [1.0]
               ]),
             o: 550.4747924804688,
             ok: false,
             p: p
           } == sn
  end

  test "match gbp4d example with first bin" do
    ldhw =
      Matrex.new([
        @it.l |> Enum.slice(-5..-1),
        @it.d |> Enum.slice(-5..-1),
        @it.h |> Enum.slice(-5..-1),
        @it.w |> Enum.slice(-5..-1)
      ])

    m =
      Matrex.new([
        @bn.l,
        @bn.d,
        @bn.h,
        @bn.w
      ])

    m1l = Matrex.column(m, 1)

    p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, m1l)

    sn = Gbp4d.gbp4d_solver_dpp(p, ldhw, m1l)
    Gbp4d.Ck.gbp4d_checkr(sn)

    assert %Gbp4d{
             bn: m1l,
             it:
               Matrex.new([
                 [-1.0, -1.0, -1.0, 0.0, -1.0],
                 [-1.0, -1.0, -1.0, 0.0, -1.0],
                 [-1.0, -1.0, -1.0, 0.0, -1.0],
                 [-1.0, -1.0, -1.0, 0.0, -1.0],
                 [2.14, 7.24, 7.24, 6.0, 4.0],
                 [3.58, 7.24, 7.24, 6.0, 4.0],
                 [4.76, 2.58, 2.58, 6.0, 4.0],
                 [243.0, 110.0, 110.0, 235.0, 258.0]
               ]),
             k:
               Matrex.new([
                 [0.0],
                 [0.0],
                 [0.0],
                 [1.0],
                 [0.0]
               ]),
             o: 216.0,
             ok: false,
             p: p
           } == sn
  end

  test "pack one item in one bin" do
    it = %{
      sku: ["16284119"],
      l: [0.1],
      d: [0.1],
      h: [0.1],
      w: [21.8]
    }

    bn = %{
      id: ["B1"],
      l: [100.0000],
      d: [100.0000],
      h: [100.0000],
      w: [22.000]
    }

    ldhw =
      Matrex.new([
        it.l,
        it.d,
        it.h,
        it.w
      ])

    m =
      Matrex.new([
        bn.l,
        bn.d,
        bn.h,
        bn.w
      ])
      |> Matrex.column(1)

    p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, m)

    sn = Gbp4d.gbp4d_solver_dpp(p, ldhw, m)

    assert %Gbp4d{
             bn: m,
             it:
               Matrex.new([
                 [0],
                 [0],
                 [0],
                 [0],
                 [0.1],
                 [0.1],
                 [0.1],
                 [21.8]
               ]),
             k:
               Matrex.new([
                 [1.0]
               ]),
             o: 0.0010000000474974513,
             ok: true,
             p: p
           } == sn

    Gbp4d.Ck.gbp4d_checkr(sn)
  end

  test "pack one item in one bin v2" do
    it = %{
      sku: ["16284119"],
      l: [0.1],
      d: [0.1],
      h: [0.1],
      w: [22.1]
    }

    bn = %{
      id: ["B1"],
      l: [100.0000],
      d: [100.0000],
      h: [100.0000],
      w: [22.000]
    }

    ldhw =
      Matrex.new([
        it.l,
        it.d,
        it.h,
        it.w
      ])

    m =
      Matrex.new([
        bn.l,
        bn.d,
        bn.h,
        bn.w
      ])
      |> Matrex.column(1)

    p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, m)

    sn = Gbp4d.gbp4d_solver_dpp(p, ldhw, m)

    assert %Gbp4d{
             bn: m,
             it:
               Matrex.new([
                 [-1],
                 [-1],
                 [-1],
                 [-1],
                 [0.1],
                 [0.1],
                 [0.1],
                 [22.1]
               ]),
             k:
               Matrex.new([
                 [0.0]
               ]),
             o: 0,
             ok: false,
             p: p
           } == sn

    Gbp4d.Ck.gbp4d_checkr(sn)
  end

  test "pack one item in one bin v3" do
    it = %{
      sku: ["16284119"],
      l: [1],
      d: [100],
      h: [100],
      w: [1]
    }

    bn = %{
      id: ["B1"],
      l: [100.0000],
      d: [100.0000],
      h: [100.0000],
      w: [22.000]
    }

    ldhw =
      Matrex.new([
        it.l,
        it.d,
        it.h,
        it.w
      ])

    m =
      Matrex.new([
        bn.l,
        bn.d,
        bn.h,
        bn.w
      ])
      |> Matrex.column(1)

    p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, m)

    sn = Gbp4d.gbp4d_solver_dpp(p, ldhw, m)

    assert %Gbp4d{
             bn: m,
             it:
               Matrex.new([
                 [0],
                 [0],
                 [0],
                 [0],
                 [1],
                 [100],
                 [100],
                 [1]
               ]),
             k:
               Matrex.new([
                 [1.0]
               ]),
             o: 1.0e4,
             ok: true,
             p: p
           } == sn

    Gbp4d.Ck.gbp4d_checkr(sn)
  end

  test "pack multiple item into one bin" do
    it = %{
      sku: ["16284119", "17850", "20136"],
      l: [40, 20, 5],
      d: [30, 10, 5],
      h: [30, 10, 5],
      w: [8, 1, 1]
    }

    bn = %{
      id: ["B1"],
      l: [100.0000],
      d: [100.0000],
      h: [100.0000],
      w: [22.000]
    }

    ldhw =
      Matrex.new([
        it.l,
        it.d,
        it.h,
        it.w
      ])

    m =
      Matrex.new([
        bn.l,
        bn.d,
        bn.h,
        bn.w
      ])
      |> Matrex.column(1)

    p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, m)

    sn = Gbp4d.gbp4d_solver_dpp(p, ldhw, m)

    assert %Gbp4d{
             bn: m,
             it:
               Matrex.new([
                 [0, 40, 60],
                 [0, 0, 0],
                 [0, 0, 0],
                 [0, 8, 9],
                 [40, 20, 5],
                 [30, 10, 5],
                 [30, 10, 5],
                 [8, 1, 1]
               ]),
             k:
               Matrex.new([
                 [1.0],
                 [1.0],
                 [1.0]
               ]),
             o: 38125.0,
             ok: true,
             p: p
           } == sn

    Gbp4d.Ck.gbp4d_checkr(sn)
  end

  test "pack multiple item into one bin v2" do
    it = %{
      sku: ["A", "A", "B", "B", "B", "B"],
      l: [60, 60, 10, 10, 10, 10],
      d: [40, 40, 7, 7, 7, 7],
      h: [40, 40, 9, 9, 9, 9],
      w: [10, 10, 0.44, 0.44, 0.44, 0.44]
    }

    bn = %{
      id: ["B1"],
      l: [100.0000],
      d: [100.0000],
      h: [100.0000],
      w: [22.000]
    }

    ldhw =
      Matrex.new([
        it.l,
        it.d,
        it.h,
        it.w
      ])

    m =
      Matrex.new([
        bn.l,
        bn.d,
        bn.h,
        bn.w
      ])
      |> Matrex.column(1)

    p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, m)

    sn = Gbp4d.gbp4d_solver_dpp(p, ldhw, m)

    Gbp4d.Ck.gbp4d_checkr(sn)

    assert %Gbp4d{
             bn: m,
             it:
               Matrex.new([
                 [0, 0, 0, 0, 0, 0],
                 [40, 0, 90, 80, 90, 80],
                 [0, 0, 9, 9, 0, 0],
                 [10, 0, 21.32, 20.88, 20.44, 20.00],
                 [60, 60, 7, 7, 7, 7],
                 [40, 40, 10, 10, 10, 10],
                 [40, 40, 9, 9, 9, 9],
                 [10, 10, 0.44, 0.44, 0.44, 0.44]
               ]),
             k:
               Matrex.new([
                 [1.0],
                 [1.0],
                 [1.0],
                 [1.0],
                 [1.0],
                 [1.0]
               ]),
             o: 194_520,
             ok: true,
             p: p
           } == sn
  end

  def gen_products(_) do
    max_weight = 100
    max_length = 200
    max_depth = 200
    max_height = 200
    step = 0.1

    Stream.flat_map(1..trunc(max_length / step), fn l ->
      Stream.flat_map(1..trunc(max_depth / step), fn d ->
        IO.inspect("l = #{l}, d = #{d}")

        Stream.flat_map(1..trunc(max_height / step), fn h ->
          Stream.map(1..trunc(max_weight / step), fn w ->
            {l, d, h, w} = {l * step, d * step, h * step, w * step}

            [
              "#{l}-#{d}-#{h}-#{w}",
              l,
              d,
              h,
              w
            ]
          end)
        end)
      end)
    end)
  end

  test "gen_product" do
    gen_products(1)
  end

  @tag timeout: :infinity
  test "crazy test" do
    bn = %{
      id: ["B1"],
      l: [100.0000],
      d: [100.0000],
      h: [100.0000],
      w: [22.000]
    }

    m =
      Matrex.new([
        bn.l,
        bn.d,
        bn.h,
        bn.w
      ])
      |> Matrex.column(1)

    Enum.each(gen_products(1), fn [sku, l, d, h, w] ->
      it = %{
        sku: [sku],
        l: [l],
        d: [d],
        h: [h],
        w: [w]
      }

      ldhw =
        Matrex.new([
          it.l,
          it.d,
          it.h,
          it.w
        ])

      try do
        p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, m)

        Gbp4d.gbp4d_solver_dpp(p, ldhw, m)
      catch
        e ->
          Logger.error(inspect("Error: #{[l, d, h, w]} #{inspect(e)}"))
      end
    end)
  end
end
