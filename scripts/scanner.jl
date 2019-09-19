using THDMMinimizer
using DelimitedFiles

function write_headers(onlya=false)
    fname_a1 = "typea1.csv"
    fname_a2 = "typea2.csv"
    fname_b = "typeb.csv"
    fname_c = "typec.csv"

    header = "m112,m122,m222,"
    header *= "lam1,lam2,lam3,lam4,lam5,"
    header *= "mu,yt,gp,g,"
    header_a = header * "nvev1,nvev2,nvev3,cbvev1,cbvev2,cbvev3\n"
    header_bc = header * "vev1,vev2,vev3\n"

    open(fname_a1, "w") do io
        write(io, header_a)
    end
    open(fname_a2, "w") do io
        write(io, header_a)
    end
    if onlya
        return
    else
        open(fname_b, "w") do io
            write(io, header_bc)
        end
        open(fname_c, "w") do io
            write(io, header_bc)
        end
    end
end

function print_current_stats(num_a1, num_a2, num_b, num_c)
    print("\u1b[1F")
    printstyled("|", color=:yellow)
    printstyled("A1 = $num_a1", color=:blue)
    printstyled("|", color=:yellow)
    printstyled("A2 = $num_a2", color=:cyan)
    printstyled("|", color=:yellow)
    printstyled("B = $num_b",   color=:green)
    printstyled("|", color=:yellow)
    printstyled("C = $num_c", color=:red)
    printstyled("|", color=:yellow)
    print("\u1b[0K")
    println()
end

function print_current_stats(num_a1, num_a2)
    print("\u1b[1F")
    printstyled("|", color=:yellow)
    printstyled("A1 = $num_a1", color=:blue)
    printstyled("|", color=:yellow)
    printstyled("A2 = $num_a2", color=:cyan)
    printstyled("|", color=:yellow)
    print("\u1b[0K")
    println()
end

function scan(;μ::Float64 = HIGGS_VEV, onlya=false)
    fname_a1 = "typea1.csv"
    fname_a2 = "typea2.csv"
    fname_b = "typeb.csv"
    fname_c = "typec.csv"

    write_headers(onlya)

    num_a1 = 0;
    num_a2 = 0;
    num_b = 0;
    num_c = 0;
    println("Starting a scan!")
    while true
        try
            nvac, cbvac, params = solve_root_eqns(μ=μ,tol=1e-5)
            vacs = [nvac, cbvac]
            vacs = [vacs; find_new_minima(params;tol=1e-5,method=Optim.BFGS())]
            type = catagorize_results(vacs)

            if type == :typea1
                deepestn = find_deepest_normal_min(vacs)
                deepestcb = find_deepest_cb_min(vacs)
                write_results_to_file(fname_a1, deepestn, deepestcb, params)
                num_a1 += 1
                if onlya
                    print_current_stats(num_a1,num_a2)
                else
                    print_current_stats(num_a1,num_a2,num_b,num_c)
                end
            elseif type == :typea2
                deepestn = find_deepest_normal_min(vacs)
                deepestcb = find_deepest_cb_min(vacs)
                write_results_to_file(fname_a2, deepestn, deepestcb, params)
                num_a2 += 1
                if onlya
                    print_current_stats(num_a1,num_a2)
                else
                    print_current_stats(num_a1,num_a2,num_b,num_c)
                end
            elseif type == :typeb && !onlya
                deepestn = find_deepest_normal_min(vacs)
                write_results_to_file(fname_b, deepestn, params)
                num_b += 1
                print_current_stats(num_a1,num_a2,num_b,num_c)
            elseif type == :typec && !onlya
                deepestcb = find_deepest_cb_min(vacs)
                write_results_to_file(fname_c, deepestcb, params)
                num_c += 1
                print_current_stats(num_a1,num_a2,num_b,num_c)
            end
        catch e
            if e isa InterruptException
                throw(e)
            else
                println("Error: $e")
                continue
            end
        end
    end
end

scan(onlya=true)
