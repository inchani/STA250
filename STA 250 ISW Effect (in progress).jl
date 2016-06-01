
path = "/Users/inchani/Desktop/UC\ Davis/My\ Courses/STA\ 250\ (AstroStatistics)/Project/"; 
map_name = "HFI_SkyMap_100_2048_R2.02_full.fits"; # Polarization map at 100 GHz
using PyCall, PyPlot
@pyimport healpy as hp
dθ = dϕ = 0.0007669903939429012;  # 5 arcmin of resolution to radian
NSIDE = 2048;
Nested = false;
I_STOKES = hp.read_map("$path$map_name", field = 0, memmap = true);
dim = length(I_STOKES);

Q_STOKES = hp.read_map("$path$map_name", field = 1, memmap = true);
U_STOKES = hp.read_map("$path$map_name", field = 2, memmap = true);

GalMapFile = "HFI_Mask_GalPlane-apo5_2048_R2.00.fits"
PtMapFile  = "HFI_Mask_PointSrc_2048_R2.00.fits"
GalMap     = hp.read_map("$path$GalMapFile", field = 3, memmap = true); # 70% sky coverage 
PtMap      = hp.read_map("$path$PtMapFile", field = 0, memmap = true);  

hp.mollview(I_STOKES*1e6, xsize = 800, title = "Intensity Map (\$\\mu\$\K\$\_{CMB}\$\)", min = -300, max = 300)

planck_map = copy(I_STOKES)
cnt = 0
for i = 1:max(length(GalMap),length(PtMap))
    if (GalMap[i] == 0.) | (PtMap[i] == 0.)
        planck_map[i] = hp.UNSEEN
        #planck_map[i] = 0.        # setting masked pixel = 0. 
        cnt += 1
    end
end

hp.mollview(planck_map*1e6, xsize = 800, title = "Intensity Map (\$\\mu\$\K\$\_{CMB}\$\)", min = -300, max = 300)

SZsource = "HFI_PCCS_SZ-union_R2.08.fits";
@pyimport healpy.fitsfunc as fitsfunc;
@pyimport healpy.pixelfunc as pixelfunc;
hdulist = fitsfunc.mrdfits("$path$SZsource",hdu=1);   #  load HDU table #1
ind = sortperm(hdulist[20][:],rev=true);              #  return indices of ojects in decreasing order of total mass

function SpheCoord2Index(θmin::Float64,θmax::Float64,ϕmin::Float64,ϕmax::Float64)
   
    #θmax > 1pi ? θmax = pi : θmax = θmax
    #θmin < 0.  ? θmin = 0. : θmin = θmin
    
    index = Array(Int64,0)
    
    if (θmax > 160pi / 180) | (θmin < 10pi / 180) 
        println("   |> Too high or low latitude... exit")
        return 0, index
    end

    
    θ = copy(θmax)
    Nθpix = 0; Nϕpix = 0; Ntot = 0
    
    Nϕpix1 = pixelfunc.ang2pix(NSIDE,θmax,ϕmax) - pixelfunc.ang2pix(NSIDE,θmax,ϕmin)
    Nϕpix2 = pixelfunc.ang2pix(NSIDE,θmin,ϕmax) - pixelfunc.ang2pix(NSIDE,θmin,ϕmin)

    i   = pixelfunc.get_all_neighbours(NSIDE,θmax,ϕmax)[6]
    l,b = pixelfunc.pix2ang(NSIDE,i)
    dϕ  = b - ϕmax
    Nϕ  = (ϕmax - ϕmin) / dϕ
    
    #println("   |> $θmin, $θmax, $ϕmin, $ϕmax, Nϕ = $Nϕ")
    println("   |> dϕ = $dϕ,  Nϕ = $Nϕ")
    if Nϕ < 1.
        println("   |> Poor resolution... exit")
        return 0, index
    end
    

    
    
    if (Nϕpix1 > 0) & (Nϕpix2 > 0)
        Nϕpix = min(Nϕpix1,Nϕpix2)
        if Nϕpix / Nϕ > 1.5 
            Nϕpix = -1
        end
    else
        Nϕpix = -1
    end

    
    if ((ϕmax > 2pi) | (ϕmin * ϕmax < 0)) 
        Nϕpix = pixelfunc.ang2pix(NSIDE,θmax,ϕmax) - pixelfunc.ang2pix(NSIDE,θmax,ϕmin)        
        N = 1; ϕ = ϕmax;    
        while (ϕ > 0.)
            i1 = pixelfunc.ang2pix(NSIDE,θmax,ϕ)
            i2 = pixelfunc.get_all_neighbours(NSIDE,θmax,ϕ)[2]
            N += 1
            l, ϕ = pixelfunc.pix2ang(NSIDE,i2)       
            if ϕ > ϕmax
                ϕ -= 2pi
            end
        end           
        N2 = copy(N-1)       
        N = 1; ϕ = 0.;
        
        while (ϕ > ϕmin)
            i1 = pixelfunc.ang2pix(NSIDE,θmax,ϕ)
            i2 = pixelfunc.get_all_neighbours(NSIDE,θmax,ϕ)[2]
            N += 1
            l, ϕ = pixelfunc.pix2ang(NSIDE,i2)
            if ϕ > 0.
                ϕ -= 2pi
            end
        end
        N1 = copy(N-1) 
        while θ >= θmin
            ϕ = ϕmin;
            ist1 = pixelfunc.ang2pix(NSIDE,θ,ϕ)
            ist2 = pixelfunc.get_all_neighbours(NSIDE,θ,ϕ)[3] # Again Pixels on North West to fully cover area.
            ind1 = Int64[i for i=ist1+1:ist1+N1]            
            ind2 = Int64[i for i=ist2+1:ist2+N1]

            index = vcat(index,ind1,ind2)                    
            ϕ = ϕmax;
            ist1 = pixelfunc.ang2pix(NSIDE,θ,ϕ)
            ist2 = pixelfunc.get_all_neighbours(NSIDE,θ,ϕ)[3] # Again Pixels on North West to fully cover area.
            ind1 = Int64[i for i=ist1-N2:ist1-1]            
            ind2 = Int64[i for i=ist2-N2:ist2-1]            
            
            index = vcat(index,ind1,ind2)        
            i = pixelfunc.get_all_neighbours(NSIDE,θ,ϕ)[4]
            θ, b = pixelfunc.pix2ang(NSIDE,i) 
            Nθpix += 2
        end
        Nϕpix = N1 + N2
        if abs( (Nϕpix-Nϕ)/Nϕ ) > 0.5
            println("   |> Found a huge discrepancy btwn Nϕpix($Nϕpix) and Nϕ... exit")
            return 0, index
        end
        println("   |> Scheme No.1 & returning ($N1+$N2=$Nϕpix, $Nθpix) array")
        return length(index), index
    end
    
    if Nϕpix < 0        
        N = 1; dind = 1;
        ϕ = ϕmin
        while (dind == 1) & (ϕ < ϕmax)
            i1 = pixelfunc.ang2pix(NSIDE,θmax,ϕ)
            i2 = pixelfunc.get_all_neighbours(NSIDE,θmax,ϕ)[6]
            dind = i2 - i1
            if dind > 0
                N += 1
            end
            l, ϕ = pixelfunc.pix2ang(NSIDE,i2)
        end           
        N1 = copy(N)
        #println("ϕ starting at $ϕ")
        N = 1; dind = 1;
        while (dind == 1) & (ϕ < ϕmax)
            i1 = pixelfunc.ang2pix(NSIDE,θmax,ϕ)
            i2 = pixelfunc.get_all_neighbours(NSIDE,θmax,ϕ)[2]
            dind = i1 - i2
            if dind > 0
                N += 1
            end
            l, ϕ = pixelfunc.pix2ang(NSIDE,i2)
        end
        #println("ϕ ending at $ϕ")        
        N2 = copy(N)
        while θ >= θmin
            ϕ = ϕmin;
            ist1 = pixelfunc.ang2pix(NSIDE,θ,ϕ)
            ist2 = pixelfunc.get_all_neighbours(NSIDE,θ,ϕ)[3] # Again Pixels on North West to fully cover area.
            ind1 = Int64[i for i=ist1+1:ist1+N1]            
            ind2 = Int64[i for i=ist2+1:ist2+N1]

            index = vcat(index,ind1,ind2)        
            #(SW, W, NW, N, NE, E, SE and S )
            
            ϕ = ϕmax;
            ist1 = pixelfunc.ang2pix(NSIDE,θ,ϕ)
            ist2 = pixelfunc.get_all_neighbours(NSIDE,θ,ϕ)[3] # Again Pixels on North West to fully cover area.
            ind1 = Int64[i for i=ist1-N2:ist1-1]            
            ind2 = Int64[i for i=ist2-N2:ist2-1]            
            
            index = vcat(index,ind1,ind2)        
            i = pixelfunc.get_all_neighbours(NSIDE,θ,ϕ)[4]
            θ, b = pixelfunc.pix2ang(NSIDE,i) 
            Nθpix += 2
        end
        Nϕpix = N1 + N2
        if abs( (Nϕpix-Nϕ)/Nϕ ) > 0.5
            println("   |> Found a huge discrepancy btwn Nϕpix($Nϕpix) and Nϕ... exit")
            return 0, index
        end        
        println("   |> Scheme No.2 & returning ($N1+$N2=$Nϕpix, $Nθpix) array")        
        return length(index), index

    end    
     
    while θ > θmin
        ϕ = ϕmin;
        ist1 = pixelfunc.ang2pix(NSIDE,θ,ϕ)
        ist2 = pixelfunc.get_all_neighbours(NSIDE,θ,ϕ)[3] # Again Pixels on North West to fully cover area.
        ind1 = Int64[i for i=ist1:ist1+Nϕpix-1]            
        ind2 = Int64[i for i=ist2:ist2+Nϕpix-1]
        
        index = vcat(index,ind1,ind2)        
        #(SW, W, NW, N, NE, E, SE and S )
        i = pixelfunc.get_all_neighbours(NSIDE,θ,ϕ)[4]
        θ, b = pixelfunc.pix2ang(NSIDE,i) 
        Nθpix += 2
    end
    println("   |> Scheme No.3 & returning ($Nϕpix, $Nθpix) array")        
    
    return length(index), index
end

#dtheta = 15. / 360 * 2pi; # 15 degrees in radian
#Nϕ, Nθ, index = SpheCoord2Index(0.3pi,0.3pi+dtheta,50dϕ,50dϕ+dtheta);
#planck_map = copy(I_STOKES)
#cnt = 0
#for i = 1:length(index)
#    planck_map[index[i]] = hp.UNSEEN
#end
#hp.mollview(planck_map*1e6, xsize = 800, title = "Intensity Map (\$\\mu\$\K\$\_{CMB}\$\)", min = -300, max = 300)

CoordSCluster = Array(Float64, length(ind), 2) # Center of Clusters in [θ, ϕ] in radian
for i = 1:length(ind)
    CoordSCluster[i,1] = pi*0.5 - hdulist[4][ind[i]]* pi / 180
    CoordSCluster[i,2] = hdulist[3][ind[i]] * pi / 180
end
dtheta = 10. / 180 * pi; # 15 degrees in radian

planck_map = copy(I_STOKES);
#for i = 1: max(length(GalMap),length(PtMap))
#    if (GalMap[i] == 0.) | (PtMap[i] == 0.)
#        planck_map[i] = 0.                  # setting masked pixel = 0. 
#    end
#end

i_sc_pixel  = Array(Int64, 0)     # a list of pixel indices of a supercluster 
id_sc       = Array(Int64, 0)     # survived index of CoordSCluster after sorting out bad samples
N_sc        = Array(Int64, 0)     # number of pixels in a region
println("---| Start clipping regions of superclusters: 20 deg x 20 deg (10 deg = $dtheta radian)")
for i = 1:60
    l, b = CoordSCluster[i,:]
    println("No. $i with angular position, (θ,ϕ) = ($l, $b)")
    N, ind = SpheCoord2Index(CoordSCluster[i,1]-dtheta,CoordSCluster[i,1]+dtheta,
    CoordSCluster[i,2]-dtheta,CoordSCluster[i,2]+dtheta)   
    if N > 0 
        planck_map[ind] = hp.UNSEEN
        id_sc           = vcat(id_sc, i)
        N_sc            = vcat(N_sc, N)
        i_sc_pixel      = vcat(i_sc_pixel, ind)
        println("   |> total N = $N")
    end
    
end

hp.mollview(planck_map*1e6, xsize = 800, title = "Intensity Map (\$\\mu\$\K\$\_{CMB}\$\)", min = -300, max = 300)

const FWHM    = 9.68 / 60. * pi / 180;   # 9.68 arcmin to radian unit (See Planck 2015 I paper)
const σ       = FWHM / √(2log(2))
const σnorm2  = 2.*σ^2.
const σlim2   = (3σ)^2.  
function GKernel(x::Float64,y::Float64,x₀::Float64,y₀::Float64)
    r2 = (x-x₀)^2. + (y-y₀)^2.
    if r2 < σlim2
        return exp( -r2 / σnorm2 )
    else
        return 0.
    end
end

const angsize    = 20.                   # width and height in degree 
const XYsize     = angsize * pi / 180;   # in radian 
Nsize            = 150;
res              = XYsize / Nsize;

function ShiftArray(X::Array{Float64,2},drow::Int64,dcol::Int64)
    Nrow  = Int(size(X)[1] / 2)
    Ncol  = Int(size(X)[2] / 2)
    rtn   = zeros(size(X)[1], size(X)[2]) 
    
    if (drow >= 0) & (dcol >= 0)
        rtn[1+drow:2Nrow,1+dcol:2Ncol] = X[1:2Nrow-drow,1:2Ncol-dcol]
        return rtn
    end
    
    if (drow <= 0) & (dcol <= 0)
        rtn[1:2Nrow+drow,1:2Ncol+dcol] = X[1-drow:2Nrow,1-dcol:2Ncol] 
        return rtn
    end
    
    if (drow >= 0) & (dcol < 0)
        rtn[1+drow:2Nrow,1:2Ncol+dcol] = X[1:2Nrow-drow,1-dcol:2Ncol] 
        return rtn
    end
    
    if (drow < 0) & (dcol >= 0)
        rtn[1:2Nrow+drow,1+dcol:2Ncol] = X[1-drow:2Nrow,1:2Ncol-dcol] 
        return rtn
    end     
end

x1d      = linspace(-XYsize*0.5,XYsize*0.5,Nsize)
y1d      = linspace(XYsize*0.5,-XYsize*0.5,Nsize)
Tmap     = Float64[ GKernel(xi,yi,0.,0.) for xi in x1d, yi in y1d];
TmapNorm =  sum(Tmap);
Tmap    /= TmapNorm;

figure(figsize=(5,5))
imshow(Tmap, origin='l', extent= [-10,10,-10,10])
colorbar();

@pyimport numpy as np
Tmin = np.min(planck_map); Tmax = np.max(planck_map)

planck_map = copy(I_STOKES) * 1e6 # Converting it to mu K scale;
Tmin = -593.5015506111085 # mu K
Tmax = 137299.59726333618 # mu K

for i = 1:max(length(GalMap),length(PtMap))
    if (GalMap[i] == 0.) | (PtMap[i] == 0.)
        planck_map[i] = hp.UNSEEN
    end
end


function pick_random_ind(index::Array{Int64,1},num::Int64)
    
    Ntot = length(index)       
    if (num == 0) | (num > Ntot)
        return index
    end

    j = 1
    ind = Array(Int64, 0)
    i = Int64( round(rand() * (Ntot-1)+1) )
    ind = vcat(ind,i)
    while j < num
        i = Int64( floor(1+ rand()*(Ntot-1)) )
        cnt = countnz(ind-i)
        #println(ind,ind-i,i)
        if cnt == j
            ind = vcat(ind,i)
            j   += 1
        end
    end

    return index[ind]
end
            
function StackImg(i0::Int64,i1::Int64,degrade=1)
    StackImage = zeros(Nsize, Nsize);   # [row, col]
    if degrade == 1
        println("degrading image is on")
    end
    println("   |> starting from $i0.")
    for i = i0:i1     # index of superclusters
        θc, ϕc = CoordSCluster[id_sc[i],:]   # Center Coord. of Supercluster
        i == 1 ? ist = 1: ist = 1 + sum(N_sc[1:i-1]) 
        percent = 20
        
        if (degrade == 1) & (Nsize*Nsize < N_sc[i])
            from = N_sc[i]
            to   = Nsize*Nsize
            println("   |> degrading img: $from to $to pixels")
            i_sc_pixel_new = pick_random_ind(i_sc_pixel[ist:ist+N_sc[i]],to)
            println("   |> degrading img: done")
            for j = 1:to
                if planck_map[i_sc_pixel_new[j]] > Tmin 
                    θ, ϕ = pixelfunc.pix2ang(NSIDE,i_sc_pixel_new[j])
                    if ϕ > ϕc + dtheta
                        ϕ -= 2pi
                    end
                    if ϕ < ϕc - dtheta
                        ϕ += 2pi
                    end
                    row_shift = Int64(round( (θ - θc) /res )) 
                    col_shift = Int64(round( (ϕ - ϕc) /res )) 
                    StackImage += ShiftArray(Tmap,row_shift, col_shift) * planck_map[i_sc_pixel_new[j]]               
                end
            end
            
            
        else     
            for j = ist:ist+N_sc[i]    # pixel index of a supercluster
                if planck_map[i_sc_pixel[j]] > Tmin # exlude the masked regions
                    θ, ϕ = pixelfunc.pix2ang(NSIDE,i_sc_pixel[j])
                    if ϕ > ϕc + dtheta
                        ϕ -= 2pi
                    end
                    if ϕ < ϕc - dtheta
                        ϕ += 2pi
                    end
                    row_shift = Int64(round( (θ - θc) /res )) 
                    col_shift = Int64(round( (ϕ - ϕc) /res )) 
                    StackImage += ShiftArray(Tmap,row_shift, col_shift) * planck_map[i_sc_pixel[j]]

                    #if Int64(round(100. * (j - ist) / N_sc[i])) > percent
                    #    println("   |> $percent % is done.")
                    #    percent += 20
                    #end
                end
            end
        end
        println("   |> No. $i is done.")
    end
    return StackImage / (i1 - i0 + 1)
end

@time StackImage = StackImg(1,49);

figure(figsize=(6,6))
imshow(StackImage, origin='l',vmin = -40, vmax = 40, extent= [-10,10,-10,10])
colorbar(); 
