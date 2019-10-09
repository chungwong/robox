defmodule Robox do
  @moduledoc """
  """

  alias Robox.{Gbp1d, Gbp4d}

  @types %{
    gbp1d: Gbp1d,
    gbp4d: Gbp4d
  }

  @spec pack(list, list, atom) :: struct
  def pack(ldhw, bn, type) do
    module = Map.get(@types, type)
    module.pack(ldhw, bn)
  end
end
