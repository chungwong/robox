defmodule Robox.Gbp4qTest do
  use ExUnit.Case
  alias Robox.{Gbp4d, Gbp4q}

  @bns [
    [06.0, 06.0, 06.0, 600.0],
    [08.0, 08.0, 06.0, 600.0],
    [09.0, 08.0, 07.0, 800.0],
    [10.0, 10.0, 10.0, 800.0],
    [22.0, 14.0, 09.0, 800.0]
  ]

  test "selects the smallest bin that fits all items" do
    sn = Gbp4q.pack([[4, 4, 4, 100], [4, 4, 4, 100]], @bns)

    assert Gbp4d.Ck.gbp4q_checkr(sn)
    assert sn.ok
    assert sn.k == [1, 1]
    # both 4x4x4 items fit the 8x8x6 bin but not the 6x6x6 one
    assert sn.f == [0, 1, 0, 0, 0]
    assert sn.o == 128.0
  end

  test "single small item selects the first bin" do
    sn = Gbp4q.pack([[1, 1, 1, 1]], @bns)

    assert Gbp4d.Ck.gbp4q_checkr(sn)
    assert sn.ok
    assert sn.f == [1, 0, 0, 0, 0]
  end

  test "no bin fits everything - picks bin maximizing fitted volume" do
    sn = Gbp4q.pack([[22, 14, 9, 100], [22, 14, 9, 100]], @bns)

    assert Gbp4d.Ck.gbp4q_checkr(sn)
    refute sn.ok
    assert sn.f == [0, 0, 0, 0, 1]
    assert Enum.sum(sn.k) == 1
  end

  test "via Robox.pack" do
    sn = Robox.pack([[4, 4, 4, 100], [4, 4, 4, 100]], @bns, :gbp4d)

    assert %Gbp4q{} = sn
    assert sn.ok
  end

  test "via Robox.pack with a single bin" do
    sn = Robox.pack([[4, 4, 4, 100]], [[6.0, 6.0, 6.0, 600.0]], :gbp4d)

    assert %Gbp4d{} = sn
    assert sn.ok
  end
end
