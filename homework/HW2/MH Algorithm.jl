
path = "/Users/inchani/Desktop/UC\ Davis/My\ Courses/STA\ 250\ (AstroStatistics)/HW/HW2"
run(`head $path/wmap_lcdm_wmap9_chains_v5/omegach2`);
nlines = countlines("$path/wmap_lcdm_wmap9_chains_v5/omegach2");

nchain = 100_000
omega_b_chain     = readdlm("$path/wmap_lcdm_wmap9_chains_v5/omegabh2")[1:nchain,2]
omega_cdm_chain   = readdlm("$path/wmap_lcdm_wmap9_chains_v5/omegach2")[1:nchain,2]
tau_reio_chain    = readdlm("$path/wmap_lcdm_wmap9_chains_v5/tau")[1:nchain,2]
theta_s_chain     = readdlm("$path/wmap_lcdm_wmap9_chains_v5/thetastar")[1:nchain,2]
A_s_109_chain     = readdlm("$path/wmap_lcdm_wmap9_chains_v5/a002")[1:nchain,2]  # <-- 10⁹ * A_s
n_s_chain         = readdlm("$path/wmap_lcdm_wmap9_chains_v5/ns002")[1:nchain,2]

full_chain        = hcat(omega_b_chain, omega_cdm_chain, tau_reio_chain, theta_s_chain, A_s_109_chain, n_s_chain);
names_chain       = [:omega_b, :omega_cdm, :tau_reio, :theta_s, :A_s_109, :n_s];

Σwmap             = cov(full_chain);

using Distributions
const lcdm_proposal_RV = MultivariateNormal(zeros(6), Σwmap)
generate_lcdm_prop(lcdm_curr)  = rand(lcdm_proposal_RV) + lcdm_curr

using HDF5
bandpowers = h5read("$path/bandpowers.h5", "bandpowers");

const σ² = (10.0/3437.75)^2;    # <---10μkarcmin noise level converted to radian pixels (1 radian = 3437.75 arcmin)
const b² = (0.0035)^2;          # <-- pixel width 0.2ᵒ ≈ 12.0 armin ≈ 0.0035 radians          

using PyCall
@pyimport pypico
const picoload = pypico.load_pico("$path/pico3_tailmonty_v34.dat")


function pico(x)     
    omega_b     = x[1]
    omega_cdm   = x[2]
    tau_reio    = x[3]
    theta_s     = x[4]
    A_s_109     = x[5]
    n_s         = x[6]
    plout::Dict{ASCIIString, Array{Float64,1}} = picoload[:get](;
        :re_optical_depth => tau_reio,
        symbol("scalar_amp(1)") =>  1e-9*A_s_109,
        :theta => theta_s,
        :ombh2 => omega_b,
        :omch2 => omega_cdm,
        symbol("scalar_spectral_index(1)") => n_s,
        :massive_neutrinos => 3.04,
        :helium_fraction => 0.25,
        :omnuh2 => 0.0,
        symbol("scalar_nrun(1)") => 0.0,
        :force     => true
    )
    clTT::Array{Float64,1} = plout["cl_TT"]
    ells   = 0:length(clTT)-1
    clTT .*= 2π ./ ells ./ (ells + 1)
    clTT[1] = 0.0
    return clTT
end

function LogPost(θ::Array{Float64,1})
    ClTT = pico(θ)
    ell  = 0:length(ClTT)-1
    Cldd = ClTT + σ² * exp(b² .* ell .* (ell + 1) ./ (8log(2)))
    rtn = 0.0
    @inbounds for l in ell[2:end]
      rtn -= log(Cldd[l+1]) * (2l+1) / 2
      rtn -= (bandpowers[l+1] / Cldd[l+1]) * (2l+1) / 2
    end
    return rtn
end


function Pprop(θ1::Array{Float64,1},θ2::Array{Float64,1})
    d = length(θ1)
    return exp(-0.5 .* transpose(θ1-θ2) * (Σwmap\ (θ1-θ2))) ./ √( (2π)^d * det(Σwmap))
end


function MH(θ₀::Array{Float64,1}, N::Int64, g::Float64)
    d      = length(θ₀)
    θcurr  = θ₀
    θ      = Array(Float64, d, N)
    αhist  = Array(Float64, N-1)    
    uhist  = Array(Float64, N-1)    
    θ[:,1] = θ₀
    Lwmap  = chol(g .* Σwmap, Val{:L}); # <------ Any problem?
    
    if N > 1 
        for i = 2:N
            θcurr = θ[:,i-1]
            #θprop = generate_lcdm_prop(θcurr)
            θprop = θcurr   .+ g .* ( Lwmap * randn(d))             
            # Assume Gaussian Proposal (Symmetric)            
            u = log(rand())
            α = min( LogPost(θprop) - LogPost(θcurr), 0.) # Log scale
            #α = min(exp(LogPost(θprop)-LogPost(θcurr)) * Pprop(θcurr,θprop) / Pprop(θprop,θcurr) ,1.)[1]            
            #θprop = θcurr   + g * rand(lcdm_proposal_RV) 
            
            
            if u <= α
                θ[:,i] = θprop[:]
            else
                θ[:,i] = θcurr[:]
            end
            αhist[i-1] = α
            uhist[i-1] = u

        end
    end
    return αhist, uhist, θ
end

θ₀ = [0.022635 , 0.113854 , 0.0889852, 0.0103905, 2.41154, 0.972309 ];
dim = 50
using PyPlot
x = linspace(1, dim, dim);

g = 0.1

αhist, uhist, θChain = MH(θ₀,dim,g);
μα = round(mean(exp(αhist))*100.) * 0.01
μu = round(mean(exp(uhist))*100.) * 0.01
omega_b_chain = transpose(θChain[1,:])
theta_s_chain = transpose(θChain[4,:])
n_s_chain     = transpose(θChain[6,:])

subplot(3,1,1)
plot(x,omega_b_chain, label = "omega_b_chain")

legend()
subplot(3,1,2)
plot(x,theta_s_chain, label = "theta_s_chain")

legend()
subplot(3,1,3)
plot(x,n_s_chain, label = "n_s_chain")
xlabel("iteration")
legend()

plot(x[2:end],exp(αhist), label = "alpha (mu=$μα), g = $g")
plot(x[2:end],exp(uhist), label = "u (mu=$μu), g = $g")
xlabel("iteration")
legend()

using PyCall, PyPlot
names_chain   = [:omega_b, :omega_cdm, :tau_reio, :theta_s, :A_s_109, :n_s]
full_chain    = transpose(θChain)
@pyimport getdist
@pyimport getdist.plots as plots
samples = getdist.MCSamples(samples=full_chain, names=names_chain)
g = plots.getSubplotPlotter()
g[:triangle_plot](samples, filled=true)
