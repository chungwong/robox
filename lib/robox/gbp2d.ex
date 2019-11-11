defmodule Robox.Gbp2d do
  alias Robox.{Gbp2d}
  defstruct [:p, :it, :bn, :k, :o, :ok]

  @gbp2d_nlvl_mkt0 2
  @gbp2d_nlvl_mkt1 3
  @gbp2d_nlvl_mkt2 5
  @gbp2d_nlvl_mkt3 8

  @moduledoc """
  inputs:
      p: vector of items profit
      ld: characterization matrix which is a 2 x N numeric matrix with each column crresponding to an item's length and depth
      m: bin’s characterization vector (m) which is a 2 x 1 numeric vector corresponding to bin’s length and depth

  w: vector of items weight
  c: weight limit constraint
  k: indicator of which items are selected
  o: weight of items selected
  ok: indicator of whether all items are selected
  """

  @doc """
  Return a column vector containing the indices of elements of X that are non-zero or satisfy a relational condition
  http://arma.sourceforge.net/docs.html#find
  """
  @spec find({%Matrex{}, atom, number}, atom, {%Matrex{}, atom, number}) :: %Matrex{}
  def find({mat1, op1, val1}, operator, {mat2, op2, val2}) do
    rst1 = find_with_op(mat1, op1, val1)
    rst2 = find_with_op(mat2, op2, val2)
    case operator do
      :and ->
        Enum.reduce(rst1, [], fn i, rst ->
          Enum.reduce(rst2, rst, fn j, rst ->
            if i != nil && j != nil && i == j do
              List.insert_at(rst, -1, [i])
            else
              rst
            end
          end)
        end)

      :or ->
        (rst1 ++ rst2)
        |> Enum.uniq()
        |> Enum.sort()
        |> Enum.map(fn i -> [i] end)
    end
    |> case do
      [] -> nil
      m -> Matrex.new(m)
    end
  end

  def find(mat, val) do
    Matrex.to_list(mat)
    |> Enum.with_index(1)
    |> Enum.flat_map(fn {m, i} ->
      if m == val do
        [i]
      else
        []
      end
    end)
  end

  def find_with_op(mat, op, val) do
    rst =
      Matrex.to_list(mat)
      |> Enum.with_index(1)

    Enum.reduce(rst, [], fn {m, i}, rt ->
      if op.(m, val) do
        List.insert_at(rt, -1, i)
      else
        rst
      end
    end)
  end

  @doc """
  Return the sub-vector of the input vector
  """
  @spec sub_vec(%Matrex{}, %Matrex{}) :: list
  def sub_vec(vec, v_indices) do
    Enum.reduce(1..v_indices[:rows], Matrex.zeros(v_indices[:rows], 1), fn i, rst ->
      Matrex.set(rst, trunc(i), 1, Enum.at(vec, i, 1))
    end)
  end

  @doc """
  Return a column vector containing the indices of elements of X that are non-zero or satisfy a relational condition
  http://arma.sourceforge.net/docs.html#find
  """
  def gbp2d_solver_dpp(p, ld, m) do
    n = p[:rows]
    q = Robox.Matrix.sort_index(p, :desc)
    bn = m
    # c++: (0,1,2), Elixir's offset starts from 1, adding 1 to start and end
    xy_ulmt = Robox.Matrix.linspace(1, 2, 2)
    # c++: (2,3,2), Elixir's offset starts from 1, adding 1 to start and end
    ld_ulmt = Robox.Matrix.linspace(3, 4, 2)

    it = Matrex.zeros(4, n)
    it = Robox.Matrix.set_rows(it, ld_ulmt, ld)
    itastr = it
    xp = Matrex.set(Matrex.zeros(4, 1), 3, 1, Matrex.at(m, 1, 1))
    xp = Matrex.set(xp, 4, 1, Matrex.at(m, 2, 1))
    nlvl = n
    # nastr: number of it oversize volume or weight limit
    nastr = gbp2d_solver_dpp_main_create_nastr(p, ld, m)
    # g: status vector: 0 open, 1 fitted, and 2 determined no fit
    g = Matrex.zeros(n, 1)
    # gastr: g asterisk: global track current best g over recursive call
    gastr = g
    v =
      Matrex.multiply(Matrex.row(it, 3), Matrex.row(it, 4))
      |> Matrex.transpose()

    vastr = min(Matrex.sum(v), Matrex.at(bn, 1, 1) * Matrex.at(bn, 2, 1))
    u = 0
    uastr = u
    k = Matrex.zeros(n, 1)
    o = 0
    n_cols = m[:cols]

    if n == 0 || n_cols == 0 do
      %Gbp2d{
        p: p,
        it: it,
        bn: bn,
        k: k,
        o: o,
        ok: true
      }
    else
      {_, it, itastr, _, gastr, _} =
        gbp2d_solver_dpp_main(bn, it, itastr, xp, q, nlvl, nastr, g, gastr, v, vastr, u, uastr)

      glmt = Robox.Matrix.find(gastr, 1)
      k = Robox.Matrix.fill(k, glmt, 1)
      o =
        Robox.Matrix.get_rows(v, glmt)
        |> Robox.Matrix.sum()

      ok = Robox.Matrix.all?(gastr, 1)
      it = Robox.Matrix.fill(it, xy_ulmt, -1)

      vals = Robox.Matrix.get_columns(itastr, glmt)
      it = Robox.Matrix.set_columns(it, glmt, vals)
      %Gbp2d{
        p: p,
        it: it,
        bn: bn,
        k: k,
        o: o,
        ok: ok
      }
    end
  end

  def gbp2d_solver_dpp_main_create_nastr(p, ld, m) do
    q = Robox.Matrix.sort_index(p)
    v = Matrex.multiply(Matrex.row(ld, 1), Matrex.row(ld, 2))
    mv = Robox.Matrix.prod(m)
    v0 = 0
    nastr = 0
    {q_size, _} = q[:size]

    {nastr, _} =
      try do
        Enum.reduce(1..q_size, {nastr, v0}, fn i, {nastr, v0} ->
          v0 = v0 + Matrex.at(v, 1, trunc(Matrex.at(q, trunc(i), 1)) + 1)

          if v0 >= mv do
            throw({:break, q_size - 1 - i, v0})
          else
            {nastr, v0}
          end

          {nastr, v0}
        end)
      catch
        {:break, nastr, v0} -> {nastr, v0}
      end
    nastr
  end

  @doc """
  create it fit via recursive fit it into bn
  """
  @spec gbp2d_solver_dpp_main(
          %Matrex{},
          %Matrex{},
          %Matrex{},
          %Matrex{},
          %Matrex{},
          number,
          number,
          %Matrex{},
          %Matrex{},
          %Matrex{},
          number,
          number,
          number
        ) ::
          {boolean, %Matrex{}, %Matrex{}, %Matrex{}, %Matrex{}, %Matrex{}}
  def gbp2d_solver_dpp_main(bn, it, itastr, xp, q, nlvl, nastr, g, gastr, v, vastr, u, uastr) do
    id = trunc(Matrex.at(q, q[:rows] - nlvl + 1, 1)) + 1
    g_size = g[:rows]

    glmt =
      Enum.reduce(1..g_size, Matrex.zeros(g_size, 1), fn i, glmt ->
        Matrex.update(glmt, i, 1, fn x -> x + 1 end)
      end)

    ok = false

    it0 = Robox.Matrix.get_columns(it, Robox.Matrix.find(g, 1))
    xp0 = xp
    kt0 = Matrex.column(it, id)

    nlmt = gbp2d_solver_dpp_main_create_nlmt(nlvl, nastr)
    ktlist = Robox.Gbp2d.It.gbp2d_it_create_ktlist(bn, it0, kt0, xp0, nlmt)
    d = Robox.Matrix.all?(gastr, 0)

    try do
      if nlvl == 1 do
        if(ktlist.n > 0) do
          it = Matrex.set_column(it, id, Matrex.column(ktlist.kt, 1))
          g = Matrex.set(g, id, 1, 1)
          u = u + Matrex.at(v, id, 1)

          {itastr, gastr, uastr} =
            if(u > uastr) do
              {it, g, u}
            else
              {itastr, gastr, uastr}
            end

          ok =
            if(uastr == vastr || Robox.Matrix.all?(gastr, 1)) do
              ok = true
              throw({:return, ok, it, itastr, xp, gastr, uastr})
            else
              ok
            end

          {ok, it, itastr, xp, gastr, uastr}
        else
          g = Matrex.set(g, id, 1, 2)
          {ok, it, itastr, xp, gastr, uastr}
        end
      else
        if(ktlist.n > 0) do
          g = Matrex.set(g, id, 1, 1)
          u = u + Matrex.at(v, id, 1)

          {ok, it, itastr, xp, gastr, uastr} =
            Enum.reduce_while(1..ktlist.n, {ok, it, itastr, xp, gastr, uastr},
              fn i, {_, it, itastr, _, gastr, uastr} ->
                it = Matrex.set_column(it, id, Matrex.column(ktlist.kt, i))
                xp = Enum.at(ktlist.xp, i - 1)

                {itastr, gastr, uastr} =
                  if u > uastr do
                    itastr = it
                    gastr = g
                    uastr = u
                    {itastr, gastr, uastr}
                  else
                    {itastr, gastr, uastr}
                  end

                if uastr == vastr || Robox.Matrix.all?(gastr, 1) do
                  ok = true
                  throw({:return, ok, it, itastr, xp, gastr, uastr})
                else
                  {ok, it, itastr, xp, gastr, uastr}
                end
                {ok, it, itastr, xp, gastr, uastr} =
                  gbp2d_solver_dpp_main(
                    bn,
                    it,
                    itastr,
                    xp,
                    q,
                    nlvl - 1,
                    nastr,
                    g,
                    gastr,
                    v,
                    vastr,
                    u,
                    uastr
                  )
                if ok == true do
                  {:halt, {ok, it, itastr, xp, gastr, uastr}}
                else
                  {:cont, {ok, it, itastr, xp, gastr, uastr}}
                end
              end)

          if !ok && nlvl < @gbp2d_nlvl_mkt2 do
            v_indices = find({gastr, &!=/2, 1}, :and, {glmt, &==/2, 1})
            rst = sub_vec(v, v_indices)
            vmiss = Matrex.sum(rst)

            if vmiss > Matrex.at(v, id, 1) do
              g = Matrex.set(g, id, 1, 2)
              u = u - Matrex.at(v, id, 1)
              it = Matrex.set(it, 1, id, 0)
              it = Matrex.set(it, 2, id, 0)
              xp = xp0
              gbp2d_solver_dpp_main(
                bn,
                it,
                itastr,
                xp,
                q,
                nlvl - 1,
                nastr,
                g,
                gastr,
                v,
                vastr,
                u,
                uastr
              )
            else
              {ok, it, itastr, xp, gastr, uastr}
            end
          else
            {ok, it, itastr, xp, gastr, uastr}
          end
        else
          g = Matrex.set(g, id, 1, 2)
          gbp2d_solver_dpp_main(
            bn,
            it,
            itastr,
            xp,
            q,
            nlvl - 1,
            nastr,
            g,
            gastr,
            v,
            vastr,
            u,
            uastr
          )
        end
      end
    catch
      {:return, ok, it, itastr, xp, gastr, uastr} ->
        {ok, it, itastr, xp, gastr, uastr}
    end
  end

  @spec gbp2d_solver_dpp_main_create_nlmt(number, number) :: number
  def gbp2d_solver_dpp_main_create_nlmt(nlvl, nastr) do
    nlmt = 1

    if nlvl < nastr do
      nlmt
    else
      cond do
        nlvl - nastr < @gbp2d_nlvl_mkt0 ->
          0

        nlvl - nastr < @gbp2d_nlvl_mkt1 ->
          5

        nlvl - nastr < @gbp2d_nlvl_mkt2 ->
          3

        nlvl - nastr < @gbp2d_nlvl_mkt3 ->
          2

        true ->
          nlmt
      end
    end
  end
end
