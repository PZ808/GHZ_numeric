# Export the radiative rho coefficient of psi4 for fixed m, summed over ell.
# Angular spectral coefficients are only an interchange representation:
# downstream seed inversion and reconstruction operate on an LGL m-mode grid.
using GeneralizedSasakiNakamura
using SpinWeightedSpheroidalHarmonics
using Printf

function main()
    a = length(ARGS)>0 ? parse(Float64,ARGS[1]) : 0.5
    lmax = length(ARGS)>1 ? parse(Int,ARGS[2]) : 8
    path = length(ARGS)>2 ? ARGS[3] : "Data/gsn_circular_a$(a)_m2.csv"
    p=10.0
    mkpath(dirname(abspath(path)))
    open(path,"w") do io
        println(io,"# M=1,mu=1,p=10,e=0,x=1; psi4rho=-Zinf; source_ell_max=$lmax")
        println(io,"# package=",pathof(GeneralizedSasakiNakamura))
        println(io,"a,r0,m,omega,L,psi4rho_real,psi4rho_imag")
        for m in (2,-2)
            summed=Dict{Int,ComplexF64}()
            w=m/(p^1.5+a)
            for l in 2:lmax
                mode=Teukolsky_pointparticle_mode(-2,l,m,0,0,a,p,0.0,1.0)
                @assert isapprox(mode.mode.omega,w; rtol=1e-13)
                sh=spin_weighted_spheroidal_harmonic(-2,l,m,a*w;method="direct")
                ls=SpinWeightedSpheroidalHarmonics.construct_all_l_in_matrix(-2,m,sh.params.N)
                for (L,b) in zip(ls,sh.coeffs/sh.normalization_const)
                    summed[L]=get(summed,L,0.0im)-mode.amplitude*b
                end
                println("a=$a l=$l m=$m Zinf=",mode.amplitude); flush(stdout)
            end
            for L in sort(collect(keys(summed)))
                c=summed[L]
                @printf(io,"%.17e,%.17e,%d,%.17e,%d,%.17e,%.17e\n",a,p,m,w,L,real(c),imag(c))
            end
            flush(io)
        end
    end
    println("Saved ",path)
end
main()
