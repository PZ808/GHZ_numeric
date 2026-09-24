# Run with an environment containing GeneralizedSasakiNakamura:
# julia --startup-file=no -O1 scripts/check_gsn_schwarzschild_mode.jl [output.csv]
# M=mu=1, a=0, p=10, e=0, x=1, ell=2, n=k=0.
# Compare complex seeds; the reference is a finite-radius fit, not an exact value.
using GeneralizedSasakiNakamura
using Printf

const root = dirname(@__DIR__)
const reference = joinpath(root, "tests/data/psi0_schwarzschild_r0_10_lmax20_mostly_minus.csv")
const output = isempty(ARGS) ? joinpath(root, "Data/gsn_schwarzschild_mode_check.csv") : ARGS[1]

function reference_modes(path)
    modes = Dict{Tuple{Int,Int},ComplexF64}()
    for line in eachline(path)
        (startswith(line, "#") || startswith(line, "ell,") || isempty(strip(line))) && continue
        c = split(strip(line), ',')
        modes[(parse(Int,c[1]),parse(Int,c[2]))] = complex(parse(Float64,c[3]),parse(Float64,c[4]))
    end
    return modes
end

function main()
    ref = reference_modes(reference)
    println("Package: ", pathof(GeneralizedSasakiNakamura))
    println("Version: ", pkgversion(GeneralizedSasakiNakamura))
    flush(stdout)
    results = Dict()
    mkpath(dirname(abspath(output)))
    open(output, "w") do io
        println(io, "# M=1,mu=1,a=0,p=10,e=0,x=1,ell=2,n=0,k=0")
        println(io, "# package_path=", pathof(GeneralizedSasakiNakamura))
        println(io, "# package_version=", pkgversion(GeneralizedSasakiNakamura))
        println(io, "method,ell,m,omega,lambda,Zinf_real,Zinf_imag")
        for method in ("isem_trapezoidal", "trapezoidal"), m in (2,-2)
            println("Computing ", method, " m=", m); flush(stdout)
            q = @time Teukolsky_pointparticle_mode(-2,2,m,0,0,0.0,10.0,0.0,1.0; method)
            results[(method,m)] = q
            @printf(io, "%s,2,%d,%.17e,%.17e,%.17e,%.17e\n",
                    method,m,q.mode.omega,q.mode.lambda,real(q.amplitude),imag(q.amplitude))
            flush(io)
        end
    end

    # In Schwarzschild both fourfold Held angular maps have eigenvalue D_l/4.
    angular = 6.0 # (ell-1)*ell*(ell+1)*(ell+2)/4 for ell=2
    for method in ("isem_trapezoidal", "trapezoidal"), m in (2,-2)
        q, reflected = results[(method,m)], results[(method,-m)]
        w = q.mode.omega
        @assert isapprox(w, m / sqrt(1000.0); rtol=1e-14)
        @assert isapprox(q.mode.lambda, 4.0; atol=1e-13)
        psi4 = -q.amplitude # rho coefficient, rho ~ -1/r
        psi4bar = (-1)^m * conj(-reflected.amplitude)
        denominator = angular^2 + 9*w^2
        seed = ((2im/w)*angular*psi4 + 6*psi4bar) / denominator
        barseed = ((2im/w)*angular*psi4bar + 6*psi4) / denominator
        expected = 2im*w^3*ref[(2,m)] / denominator
        expected_bar = (-1)^m * conj(2im*(-w)^3*ref[(2,-m)] / denominator)
        error = abs(seed-expected)/abs(expected)
        barerror = abs(barseed-expected_bar)/abs(expected_bar)
        inferred_psi0 = denominator*seed/(2im*w^3)
        @printf("%s m=%d Zinf=%.16e %+.16ei\n",method,m,real(q.amplitude),imag(q.amplitude))
        println("  f from psi4 = ",seed,"; reference = ",expected)
        println("  inferred psi0 = ",inferred_psi0)
        @printf("  relative f error=%.10e; fbar error=%.10e\n",error,barerror)
        # Regression tolerance for this particular finite-radius reference fit.
        @assert max(error,barerror) < 2e-5 "Complex seed disagrees with Schwarzschild reference"
        # Check the original coupled equation, not just its eliminated form.
        @assert abs(angular*seed+3im*w*barseed-(2im/w)*psi4) < 1e-13
    end
    for m in (2,-2)
        a = results[("isem_trapezoidal",m)].amplitude
        b = results[("trapezoidal",m)].amplitude
        err = abs(a-b)/abs(a)
        @printf("ISEM versus legacy amplitude m=%d relative error=%.10e\n",m,err)
        @assert err < 1e-7 "Independent radial paths disagree"
    end
    println("PASS. Raw amplitudes saved to ",output)
end

main()
