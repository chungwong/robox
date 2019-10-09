defmodule Robox.Gbp4d do
  alias Robox.{Gbp1d, Gbp4d, GbpU, Matrix}
  defstruct [:p, :it, :bn, :k, :o, :ok]

  # gbp4d global setup on parameter
  # when <  2 items search over all in ktlist w.r.t orientation + extreme point
  @gbp4d_nlvl_mkt0 2
  # when <  3 items search best   5 in ktlist w.r.t orientation + extreme point
  @gbp4d_nlvl_mkt1 3
  # when <  5 items search best   3 in ktlist w.r.t orientation + extreme point
  @gbp4d_nlvl_mkt2 5
  # when <  8 items search best   2 in ktlist w.r.t orientation + extreme point
  @gbp4d_nlvl_mkt3 8
  # @gbp4d_nlvl_mkt4 = ; # when >= 8 items search best   1 in ktlist w.r.t orientation + extreme point

  @spec pack(%Matrex{}, %Matrex{}) :: %Gbp4d{}
  def pack(%Matrex{} = ldhw, %Matrex{} = bn) do
    if bn[:cols] == 1 do
      p = Gbp4d.gbp4d_solver_dpp_prep_create_p(ldhw, bn)
      gbp4d_solver_dpp(p, ldhw, bn)
    end
  end

  @spec pack(list, list) :: %Gbp4d{}
  def pack(ldhw, bn) do
    ldhw =
      Matrex.new(ldhw)
      |> Matrex.transpose()

    bn =
      Matrex.new(bn)
      |> Matrex.transpose()

    pack(ldhw, bn)
  end

  @doc """
  solve gbp4d via extreme point heuristic and best information score fit strategy
  """
  @spec gbp4d_solver_dpp(%Matrex{}, %Matrex{}, %Matrex{}) :: %Gbp4d{}
  def gbp4d_solver_dpp(p, ldhw, m) do
    # n: number of it
    {n, _} = p[:size]

    # q: it fit sequence
    q = Matrix.sort_index(p, :desc)

    bn = m
    it = Matrex.zeros(8, n)

    # it - x, y, z, w on row 0, 1, 2, 3
    # one-based index
    xyzw_ulmt = Matrix.linspace(1, 4, 4)
    ld_ulmt = Matrix.linspace(5, 8, 4)

    it = Matrix.set_rows(it, ld_ulmt, ldhw)

    # itastr: it asterisk: global track current best it over recursive call
    itastr = it

    xp =
      Matrex.zeros(8, 1)
      |> Matrex.set(5, 1, bn[1])
      |> Matrex.set(6, 1, bn[2])
      |> Matrex.set(7, 1, bn[3])
      |> Matrex.set(8, 1, bn[4])

    # nlvl: number of it open to fit
    nlvl = n

    # # nastr: number of it oversize volume or weight limit
    nastr = gbp4d_solver_dpp_main_create_nastr(p, ldhw, m)

    g = Matrex.zeros(n, 1)

    # # gastr: g asterisk: global track current best g over recursive call
    gastr = g

    # # v: volume of it - l d h in objective and constraint while w in constraint only
    v =
      Matrex.multiply(Matrex.row(it, 5), Matrex.row(it, 6))
      |> Matrex.multiply(Matrex.row(it, 7))
      |> Matrex.transpose()

    # # vastr: v asterisk: maximum volume can be achieved via fit
    vastr = min(Matrix.sum(v), bn[1] * bn[2] * bn[3])

    # # u: sum of v volume fitted into bn
    u = 0

    # # uastr: u asterisk: global track current best u over recursive call
    uastr = u

    # # k: selection indicator <vector>
    k = Matrex.zeros(n, 1)

    # # o: objective achievement <double>
    o = 0

    # main fit corner case when n == 0
    if n == 0 || m == nil || m[:size] == {0, 0} do
      %Gbp4d{
        p: p,
        it: it,
        bn: bn,
        k: k,
        o: o,
        ok: true
      }
    else
      # main fit recursive call - should create a class gbp4d_allinfo hold all info?
      {_ok, it, itastr, _xp, gastr, _uastr} =
        gbp4d_solver_dpp_main(bn, it, itastr, xp, q, nlvl, nastr, g, gastr, v, vastr, u, uastr)

      # ok via gbp4d_solver_dpp(): can all it fit into bn?
      # ok via gbp4d_solver_dpp_main() hold a different meaning for recursive purpose
      # ok via gbp4d_solver_dpp_main(): can all it fit into bn or a subset of it make full use of bn 100% volume utilization?

      glmt = Matrix.find(gastr, 1)
      k = Matrix.fill(k, glmt, 1)

      o =
        Matrix.get_rows(v, glmt)
        |> Matrix.sum()

      ok = Matrix.all?(gastr, 1)

      # # # it via itastr and gastr: flag no fit with (x, y) = (-1, -1, -1) instead of (0, 0, 0)
      it = Matrix.fill(it, xyzw_ulmt, -1)

      vals = Matrix.get_columns(itastr, glmt)
      it = Matrix.set_columns(it, glmt, vals)

      %Gbp4d{
        p: p,
        it: it,
        bn: bn,
        k: k,
        o: o,
        ok: ok
      }
    end
  end

  # @spec gbp4d_solver_dpp_main(%Matrex{}, %Matrex{}, %Matrex{}, %Matrex{}, %Matrex{}, number, number, %Matrex{}, %Matrex{}, %Matrex{}, number, number, number) :: {boolean, %Matrex{}}
  def gbp4d_solver_dpp_main(
        bn,
        it,
        itastr,
        xp,
        q,
        nlvl,
        nastr,
        g,
        gastr,
        v,
        vastr,
        u,
        uastr,
        flag \\ false
      ) do
    # index of it to fit
    id = trunc(q[q[:rows] - nlvl + 1]) + 1

    glmt =
      Matrex.zeros(g[:rows], 1)
      |> Matrix.get_rows(Matrix.find(g, 0))
      |> Matrex.add(1)

    ok = false

    # create ktlist
    it0 = Matrix.get_columns(it, Matrix.find(g, 1))
    xp0 = xp
    kt0 = Matrex.column(it, id)
    nlmt = gbp4d_solver_dpp_main_create_nlmt(nlvl, nastr)

    ktlist = Gbp4d.It.gbp4d_it_create_ktlist(bn, it0, xp0, kt0, nlmt)

    try do
      # main
      if nlvl == 1 do
        # fit final it in queue
        if ktlist.n > 0 do
          it = Matrex.set_column(it, id, Matrex.column(ktlist.kt, 1))
          g = Matrex.set(g, id, 1, 1)
          u = u + v[id]

          {itastr, gastr, uastr} =
            if u > uastr do
              {it, g, u}
            else
              {itastr, gastr, uastr}
            end

          ok =
            if uastr == vastr || Matrix.all?(gastr, 1) do
              true
            else
              ok
            end

          {ok, it, itastr, xp, gastr, uastr}
        else
          g = Matrex.set(g, id, 1, 2)
          {ok, it, itastr, xp, gastr, uastr}
        end
      else
        # do recursive call

        if ktlist.n > 0 do
          g = Matrex.set(g, id, 1, 1)

          u = u + v[id]

          {ok, it, itastr, xp, gastr, uastr} =
            Enum.reduce_while(1..ktlist.n, {ok, it, itastr, xp, gastr, uastr}, fn i,
                                                                                  {_, it, itastr,
                                                                                   _, gastr,
                                                                                   uastr} ->
              it = Matrex.set_column(it, id, Matrex.column(ktlist.kt, i))
              xp = Enum.at(ktlist.xp, i - 1)

              {itastr, gastr, uastr} =
                if u > uastr do
                  {it, g, u}
                else
                  {itastr, gastr, uastr}
                end

              if uastr == vastr || Matrix.all?(gastr, 1) do
                ok = true
                throw({:return, {ok, it, itastr, xp, gastr, uastr}})
              else
                {ok, it, itastr, xp, gastr, uastr} =
                  gbp4d_solver_dpp_main(
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
                    uastr,
                    true
                  )

                if ok do
                  {:halt, {ok, it, itastr, xp, gastr, uastr}}
                else
                  {:cont, {ok, it, itastr, xp, gastr, uastr}}
                end
              end
            end)

          # what if skip this it? can later it combined be better?
          if !ok && nlvl < @gbp4d_nlvl_mkt2 do
            vmiss =
              case Matrix.find({gastr, &!=/2, 1}, :and, {glmt, &==/2, 1}) do
                nil ->
                  0.0

                r ->
                  Matrix.get_rows(v, r)
                  |> Matrix.sum()
              end

            if vmiss > v[id] do
              g = Matrex.set(g, id, 1, 2)
              u = u - v[id]

              it =
                it
                |> Matrex.set(1, id, 0)
                |> Matrex.set(2, id, 0)
                |> Matrex.set(3, id, 0)
                |> Matrex.set(4, id, 0)

              xp = xp0

              gbp4d_solver_dpp_main(
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

          gbp4d_solver_dpp_main(
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
      {:return, {ok, it, itastr, xp, gastr, uastr}} ->
        {ok, it, itastr, xp, gastr, uastr}
    end
  end

  @spec gbp4d_solver_dpp_main_create_nastr(%Matrex{}, %Matrex{}, %Matrex{}) :: number
  def gbp4d_solver_dpp_main_create_nastr(p, ldhw, m) do
    # init
    q = Matrix.sort_index(p, :desc)

    v =
      Matrex.multiply(Matrex.row(ldhw, 1), Matrex.row(ldhw, 2))
      |> Matrex.multiply(Matrex.row(ldhw, 3))
      |> Matrex.transpose()

    w = Matrex.transpose(Matrex.row(ldhw, 4))

    mv = Matrix.prod(Matrix.subvec(m, 1, 3))

    mw = m[4]

    # main
    v0 = 0.00
    w0 = 0.00
    nastr = 0

    {_, _, nastr} =
      Enum.reduce_while(1..q[:rows], {v0, w0, nastr}, fn i, {v0, w0, nastr} ->
        qi = trunc(q[i] + 1)
        v0 = v0 + v[qi]
        w0 = w0 + w[qi]

        if v0 >= mv || w0 >= mw do
          nastr = q[:rows] - i
          {:halt, {v0, w0, nastr}}
        else
          {:cont, {v0, w0, nastr}}
        end
      end)

    nastr
  end

  def gbp4d_solver_dpp_prep_create_p(ldhw, m) do
    # init
    p = Matrex.zeros(ldhw[:cols], 1)

    if ldhw == nil || ldhw[:cols] == 0 || m == nil do
      p
    else
      # main
      q = Matrex.zeros(3, ldhw[:cols])

      # create it w cluster w.r.t gbp1d
      v =
        Matrex.multiply(Matrex.row(ldhw, 1), Matrex.row(ldhw, 2))
        |> Matrex.multiply(Matrex.row(ldhw, 3))

      w =
        Matrex.divide(Matrex.row(ldhw, 4), 0.25)
        |> Matrex.apply(:ceil)

      c = Float.floor(m[4] / 0.25)

      wfit = Gbp1d.gbp1d_solver_dpp(v, w, c)

      q =
        Enum.reduce(1..ldhw[:cols], q, fn i, q ->
          if wfit.k[i] == 1 do
            Matrex.set(q, 1, i, 1.0)
          else
            Matrex.set(q, 1, i, 0.0)
          end
        end)

      # create it max(l, d, h) cluster w.r.t bn max(l, d, h) into per 0.25 inch
      mldh =
        Matrix.subvec(m, 1, 3)
        |> Matrex.max()

      ildh =
        Matrix.get_rows(ldhw, [1, 2, 3])
        |> Matrix.max()

      nldh = Matrix.linspace(0, mldh, Float.ceil(mldh * 4.0))

      q =
        Enum.reduce(1..ldhw[:cols], q, fn i, q ->
          Matrex.set(
            q,
            2,
            i,
            Matrix.find(nldh, &<=/2, ildh[i])
            |> Matrex.subtract(1)
            |> Matrix.sum()
          )
        end)

      row =
        Matrex.multiply(Matrex.row(ldhw, 1), Matrex.row(ldhw, 2))
        |> Matrex.multiply(Matrex.row(ldhw, 3))
        |> Matrex.divide(ildh)

      q = Matrix.set_rows(q, Matrex.new([[3]]), row)

      # sort index after sort index via rows so that it with higher value of p can get fit into bn earlier in queue
      GbpU.sort_index_via_rows(q, Matrix.linspace(0, 2, 3))
      |> Matrix.sort_index()
    end
  end

  def gbp4d_solver_dpp_main_create_nlmt(nlvl, nastr) do
    nlmt = 1

    if nlvl <= nastr do
      nlmt
    else
      cond do
        nlvl - nastr < @gbp4d_nlvl_mkt0 ->
          0

        nlvl - nastr < @gbp4d_nlvl_mkt1 ->
          5

        nlvl - nastr < @gbp4d_nlvl_mkt2 ->
          3

        nlvl - nastr < @gbp4d_nlvl_mkt3 ->
          2

        true ->
          nlmt
      end
    end
  end
end
