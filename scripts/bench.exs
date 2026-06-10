# Timing benchmark across item counts.
# Usage: mix run scripts/bench.exs

:rand.seed(:exsss, {7, 49, 91})

rand_dim = fn max -> Float.round(0.5 + :rand.uniform() * max, 2) end

bench = fn n ->
  items =
    Enum.map(1..n, fn _ ->
      [rand_dim.(20.0), rand_dim.(20.0), rand_dim.(20.0), Float.round(0.5 + :rand.uniform() * 10.0, 2)]
    end)

  bin = [40.0, 35.0, 30.0, 500.0]

  {us, sn} = :timer.tc(fn -> Robox.pack(items, [bin], :gbp4d) end)

  fitted = Enum.sum(sn.k)
  IO.puts("n=#{n}: #{Float.round(us / 1000, 1)} ms, fitted #{fitted}/#{n}, ok=#{sn.ok}")
end

Enum.each([5, 10, 15, 20, 30, 50, 80], bench)

# the production crash case shape: repeated identical items
items = List.duplicate([7.24, 7.24, 2.58, 110.0], 20)
{us, sn} = :timer.tc(fn -> Robox.pack(items, [[22.0, 14.0, 9.0, 800.0]], :gbp4d) end)
IO.puts("20 identical items: #{Float.round(us / 1000, 1)} ms, fitted #{Enum.sum(sn.k)}/20, ok=#{sn.ok}")
