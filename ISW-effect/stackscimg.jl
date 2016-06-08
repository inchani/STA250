include("lib.jl")

Coord_i   = 1
CoordName = ["Coord_SZ","Coord_GR08_sc","Coord_GR08_void"]
Nsample = 1
discrepancy = 1.
const dtheta_deg = 15.


using PyCall, PyPlot
@pyimport healpy as hp
@pyimport healpy.fitsfunc as fitsfunc;
@pyimport healpy.pixelfunc as pixelfunc;
@pyimport numpy as np;

println("---| Reading a polarization map at 100 GHz")
path = "/home/inchani/STA250/";
map_name = "HFI_SkyMap_100_2048_R2.02_full.fits"; # Polarization map at 100 GHz

dθ = dϕ = 0.0007669903939429012;  # 5 arcmin of resolution in radian
NSIDE = 2048;
Nested = false;
I_STOKES = hp.read_map("$path$map_name", field = 0, memmap = true);
dim = length(I_STOKES);


const dtheta = dtheta_deg / 180 * pi;   # 15 degrees in radian
const FWHM    = 30. / 60. * pi / 180;   # 30. arcmin to radian unit (See Planck 2013 ISW paper)
const σ       = FWHM / 2.355            # 2.355 ~ 2√(2log(2))
const σnorm2  = 2.*σ^2.
const σlim2   = (3σ)^2.  

const angsize    = copy(dtheta_deg*2)        # width and height in degree 
const XYsize     = angsize * pi / 180;       # in radian 
const res        = 6. / 60. * pi / 180;  
        # 6 arcmin of pixel size in radian unit (See Planck 2013 ISW paper)
const Nsize      = Int64(round(XYsize/res)); # 300

const Tmin = -593.5015506111085; # mu K  
const Tmax =  709.0113358572125; # mu K 
# min and max temperature if we take out 40 % of sky as foreground


#In a galactic mask map, I ruled out " mask value = 0" pixels in 70% coverage case
println("---| Reading mask maps")
GalMapFile = "HFI_Mask_GalPlane-apo5_2048_R2.00.fits"
PtMapFile  = "HFI_Mask_PointSrc_2048_R2.00.fits"
GalMap     = hp.read_map("$path$GalMapFile", field = 3, memmap = true); # 70% sky coverage 
PtMap      = hp.read_map("$path$PtMapFile", field = 0, memmap = true);  


println("---| Reading coordinates of "CoordName[Coord_i])

#CoordSZ       = np.load("$path""Coord_SZ.npy");
#CoordGR08SC   = np.load("$path""Coord_GR08_sc.npy");
#CoordGR08Void = np.load("$path""Coord_GR08_void.npy");

#####################################
CoordSCluster = np.load("$path"CoordName[Coord_i]".npy");
#####################################

i_sc_pixel  = Array(Int64, 0)     # a list of pixel indices of a supercluster 
id_sc       = Array(Int64, 0)     # survived index of CoordSCluster after sorting out bad samples
N_sc        = Array(Int64, 0)     # number of pixels in a region

println("---| Start clipping regions of superclusters: 30 deg x 30 deg (15 deg = $dtheta radian)")
i = 0; ntry = 1;

if Nsample > length(CoordSCluster[:,1])
    Nsample = length(CoordSCluster[:,1])
end
while (i < Nsample) & (ntry < length(CoordSCluster[:,1]))
    l, b = CoordSCluster[ntry,:]
    θ = l * 180 / pi
    ϕ = b * 180 / pi
    #if (( θ > 15. ) & (θ < 75.)) | (( θ > 115. ) & (θ < 165.))
    N, ind = SpheCoord2Index(CoordSCluster[ntry,1]-dtheta,CoordSCluster[ntry,1]+dtheta,
    CoordSCluster[ntry,2]-dtheta,CoordSCluster[ntry,2]+dtheta;discrepancy=discrepancy)  
    if N > 0 
        i +=1        
        println("   |>> No. $i out of $ntry with angular position, (θ,ϕ) = ($θ, $ϕ)")                
        id_sc           = vcat(id_sc, i)
        N_sc            = vcat(N_sc, N)
        i_sc_pixel      = vcat(i_sc_pixel, ind)
        println("   |>> total N = $N")
    end
    #end
    ntry +=1
end

Nstacked = copy(i)


println("---| Constructing a Gaussian kernel of $Nsize by $Nsize map")
x1d      = linspace(-XYsize*0.5,XYsize*0.5,Nsize)
y1d      = linspace(XYsize*0.5,-XYsize*0.5,Nsize)
Tmap     = Float64[ GKernel(xi,yi,0.,0.) for xi in x1d, yi in y1d];
Umap     = Float64[ GKernel(xi,yi,0.,0.) > 0.? 1.: 0. for xi in x1d, yi in y1d];
TmapNorm =  sum(Tmap);
Tmap    /= TmapNorm;



# masking galactic foregrounds and point sources
planck_map = copy(I_STOKES) * 1e6 # Converting it to mu K scale;


for i = 1:max(length(GalMap),length(PtMap))
    if (GalMap[i] == 0.) | (PtMap[i] == 0.)
        planck_map[i] = hp.UNSEEN
    end
end


name = CoordName[Coord_i]

println("---| Start Stacking images")
@time StackImage = StackImg(1,Nsample; degrade = 0);
np.save("$path""stacked$Nstacked""img_$name",StackImage)
