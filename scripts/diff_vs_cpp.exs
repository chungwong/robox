# Differential test: random cases through the original C++ gbp4d solver and
# this port, comparing p / k / o / ok.
#
# Usage: MAX_ITEMS=16 mix run scripts/diff_vs_cpp.exs [ncases] [seed]
#
# Requires the compiled reference driver at /tmp/gbp_ref_src/diff_driver,
# built from the original C++ sources (see scripts/diff_driver.cpp):
#
#   mkdir -p /tmp/gbp_ref_src && cd /tmp/gbp_ref_src
#   cp <gbp checkout>/src/{gbp4d,gbp4d_xp,gbp4d_it,gbp1d,gbp_u}.{cpp,h} .
#   sed -i 's|#include <RcppArmadillo.h>|#include <armadillo>|; s|using namespace Rcpp;||' *.cpp *.h
#   cp <robox checkout>/scripts/diff_driver.cpp .
#   g++ -O2 -o diff_driver diff_driver.cpp gbp4d.cpp gbp4d_xp.cpp gbp4d_it.cpp gbp1d.cpp gbp_u.cpp -larmadillo

driver = "/tmp/gbp_ref_src/diff_driver"

{ncases, seed} =
  case System.argv() do
    [n, s] -> {String.to_integer(n), String.to_integer(s)}
    [n] -> {String.to_integer(n), 42}
    _ -> {200, 42}
  end

:rand.seed(:exsss, {seed, seed * 7, seed * 13})

rand_dim = fn max -> Float.round(0.5 + :rand.uniform() * max, 2) end

cases =
  Enum.map(1..ncases, fn _ ->
    n = :rand.uniform(String.to_integer(System.get_env("MAX_ITEMS", "8")))

    # mix of small and large bins; weight cap sometimes binding
    bin = [
      rand_dim.(40.0),
      rand_dim.(40.0),
      rand_dim.(40.0),
      Float.round(5.0 + :rand.uniform() * 50.0, 2)
    ]

    items =
      Enum.map(1..n, fn _ ->
        [
          rand_dim.(20.0),
          rand_dim.(20.0),
          rand_dim.(20.0),
          Float.round(0.5 + :rand.uniform() * 20.0, 2)
        ]
      end)

    {items, bin}
  end)

input =
  [
    Integer.to_string(ncases)
    | Enum.map(cases, fn {items, [bl, bd, bh, bw]} ->
        n = length(items)

        [
          "#{n} #{bl} #{bd} #{bh} #{bw}"
          | Enum.map(items, fn [l, d, h, w] -> "#{l} #{d} #{h} #{w}" end)
        ]
        |> Enum.join("\n")
      end)
  ]
  |> Enum.join("\n")

File.write!("/tmp/diff_cases.txt", input <> "\n")

{out, 0} = System.shell("#{driver} < /tmp/diff_cases.txt")

ref =
  out
  |> String.split("case ", trim: true)
  |> Enum.map(fn block ->
    lines = String.split(block, "\n", trim: true)
    [_idx | rest] = lines

    Enum.reduce(rest, %{}, fn line, acc ->
      case String.split(line, ":", parts: 2) do
        ["p", v] -> Map.put(acc, :p, v |> String.split() |> Enum.map(&elem(Float.parse(&1), 0)))
        ["k", v] -> Map.put(acc, :k, v |> String.split() |> Enum.map(&String.to_integer/1))
        ["o", v] -> Map.put(acc, :o, v |> String.trim() |> Float.parse() |> elem(0))
        ["ok", v] -> Map.put(acc, :ok, String.trim(v) == "1")
        _ -> acc
      end
    end)
  end)

stats = %{p: 0, k: 0, o: 0, ok: 0, checkr: 0}

stats =
  Enum.zip(cases, ref)
  |> Enum.with_index()
  |> Enum.reduce(stats, fn {{{items, bin}, ref}, idx}, stats ->
    sn = Robox.pack(items, [bin], :gbp4d)

    p_ok = Enum.map(sn.p, &(&1 / 1)) == ref.p
    k_ok = sn.k == ref.k
    o_ok = abs(sn.o - ref.o) < 1.0e-6 * max(1.0, abs(ref.o))
    okok = sn.ok == ref.ok
    ck = Robox.Gbp4d.Ck.gbp4d_checkr(sn)

    unless p_ok and k_ok and o_ok and okok and ck do
      IO.puts(
        "case #{idx}: p=#{p_ok} k=#{k_ok} o=#{o_ok} ok=#{okok} checkr=#{ck}\n" <>
          "  items=#{inspect(items)}\n  bin=#{inspect(bin)}\n" <>
          "  ref: k=#{inspect(ref.k)} o=#{ref.o} ok=#{ref.ok}\n" <>
          "  got: k=#{inspect(sn.k)} o=#{sn.o} ok=#{sn.ok}"
      )
    end

    %{
      p: stats.p + if(p_ok, do: 0, else: 1),
      k: stats.k + if(k_ok, do: 0, else: 1),
      o: stats.o + if(o_ok, do: 0, else: 1),
      ok: stats.ok + if(okok, do: 0, else: 1),
      checkr: stats.checkr + if(ck, do: 0, else: 1)
    }
  end)

IO.puts("\n#{ncases} cases, mismatches: #{inspect(stats)}")
