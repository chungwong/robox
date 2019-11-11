defmodule Robox.MixProject do
  use Mix.Project

  def project do
    [
      app: :robox,
      version: "0.1.0",
      elixir: "~> 1.6",
      package: package(),
      elixirc_paths: elixirc_paths(Mix.env()),
      start_permanent: Mix.env() == :prod,
      deps: deps()
    ]
  end

  defp elixirc_paths(:test), do: ["lib", "test/support"]
  defp elixirc_paths(_), do: ["lib"]

  defp package do
    [
      description: "Bin Packing for Elixir",
      files: [
        "lib",
        "mix.exs",
        "README.md",
        "CHANGELOG.md",
        ".formatter.exs"
      ],
      maintainers: [
        "Chung Wong"
      ],
      licenses: ["MIT"],
      links: %{}
    ]
  end

  # Run "mix help compile.app" to learn about applications.
  def application do
    [
      extra_applications: [:logger]
    ]
  end

  # Run "mix help deps" to learn about dependencies.
  defp deps do
    [
      {:matrex, "~> 0.6.0"}
    ]
  end
end
