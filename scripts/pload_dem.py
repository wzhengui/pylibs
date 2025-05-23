#!/usr/bin/env python3
'''
load bathymetry for NWM model grid 
Remove mpi4py if you only use serial
'''
from pylib import *
from mpi4py import MPI
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
grd='hgrid.ll'  #grid name (.ll or .gr3)
grdout='hgrid.ll.new' #grid name with depth loaded

#DEM informat
sdir=r'./DEM'    #directory of DEM data
format_dem='tif' #format of DEM data (npz, tif, tiff); positions are not needed for *.npz files 
reverse_sign=1   #invert depth sign
datum_shift =0   #change vertical datum (e.g. -0.258 for navd to ngvd (msl) in ChesBay)

#parameter
headers=('Bayonne','New_Arthur','CT_River','NY_TACC','Hudson_River','Long_Island','Raritan_Bay_River','MA_TACC','Toms_River')
positions=(0,0,0,0,0,0,0,0,0)  #0: cell center;  1: cell corner for DEM file (Property->AREA_OR_POINT=Point to find out this)
#headers=("etopo1","crm_3arcs","cdem13_","dem_continetalus_southcarolina","North_Carolina_USGS_3m",
#         "al_ll","nc_ll","fl_ll","gulf_1_dem_usgs","gulf_3_demcombined_ll","ge_ll","sc_ll",
#         "cb_ll","db_ll","new_england_topobathy_dem_3m_dd","Tile3_R3_DEMv2","cb_bay_dem_v3.1_ll") #FOR STOFS3D
#regions=("min_5m_ll.reg","SabinePass.reg","BergenPoint.reg","Washington_3.reg",
#         "Elk_river.reg","Hudson_river.reg","James_river.reg","NorthEast_river.reg",
#         "Rappahannock_river.reg","Susquehanna_river.reg","York_river.reg",
#         "Androscoggin_Kennebec_rivers.reg","Merrimack_river.reg","Patuxent_river.reg",
#         "Penobscot_river.reg","Saco_river.reg","StCroix_river.reg") #regions for modifying depth
#rvalues=(5,7,5,15,2,16,14,5,6,10,10,3,3,5,5,3,5) #minimum depth in regions (note: region will be skipped if not exist)

#resource requst 
ibatch=0     #0: serial mode;   1: parallel mode (for serial node, walltime/nnode/ppn are optional)
walltime='00:10:00'; nnode=1;  ppn=4
#hpc: femto, hurricane, bora, vortex, potomac, james, frontera, levante, stampede2
#ppn:   32,       8,     8,    12,       12,     20,     56,      128,      48

#optional: (frontera,levante,stampede2)
qname   ='compute'         #partition name
account ='TG-OCE140024'    #stampede2: NOAA_CSDL_NWI,TG-OCE140024; levante: gg0028 
qnode   =None              #specify node name, or default qnode based on HOST will be used

jname='load_dem'; scrout='screen.out'; bdir=os.path.abspath(os.path.curdir)
#-----------------------------------------------------------------------------
#on front node: 1). submit jobs first (qsub), 2) running parallel jobs (mpirun) 
#-----------------------------------------------------------------------------
if ibatch==0: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run locally
if os.getenv('job_on_node')==None:
   if os.getenv('param')==None: fmt=0; bcode=sys.argv[0]; os.environ['qnode']=get_qnode(qnode)
   if os.getenv('param')!=None: fmt=1; bdir,bcode=os.getenv('param').split(); os.chdir(bdir)
   scode=get_hpc_command(bcode,bdir,jname,qnode,nnode,ppn,walltime,scrout,fmt=fmt,qname=qname,account=account)
   print(scode); os.system(scode); os._exit(0)

#-----------------------------------------------------------------------------
#on computation node
#-----------------------------------------------------------------------------
bdir=os.getenv('bdir'); os.chdir(bdir) #enter working dir
if ibatch==0: nproc=1; myrank=0
if ibatch==1: comm=MPI.COMM_WORLD; nproc=comm.Get_size(); myrank=comm.Get_rank()
if myrank==0: t0=time.time()

if 'regions' not in locals(): regions=None; rvalues=None
if format_dem!='npz' and len(headers)!=len(positions): sys.exit('different size: headers, positions') 
#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
#get all DEM files and distribute jobs
fnames0=array([i for i in os.listdir(sdir) if i.endswith('.'+format_dem)])

#filter with headers, and sort by id numbers
fnames_sort=[]; ps0=[]
for header,position in zip(headers,positions):
    fnames_sub=array([i for i in fnames0 if i.startswith(header)]); psi='center' if position==0 else 'corner'
    if len(fnames_sub)==1: fnames_sort.extend(fnames_sub); ps0.append(psi); continue

    #get id number 
    if format_dem=='npz':
       fid=array([i.replace('tif.','').replace('.','_')[len(header):].split('_')[-2] for i in fnames_sub]).astype('int')
    elif format_dem=='tif':
       fid=[i for i in fnames_sub] #can add order number in the DEM name to sort it
    fnames_sort.extend(fnames_sub[argsort(fid)]); ps0.extend(tile(psi,len(fid)))
fnames_sort=array(fnames_sort); ps0=array(ps0)

#distribute jobs
fnames=[]; inum=[]; ps=[]
for m,[fname,psi] in enumerate(zip(fnames_sort,ps0)):
    if m%nproc==myrank: fnames.append(fname); inum.append(m); ps.append(psi)
fnames=array(fnames); ps=array(ps)

#read hgrid 
if grd.endswith('npz'):
   gd=loadz(grd).hgrid
elif grd.endswith('gr3') or grd.endswith('ll'):
   gd=read_schism_hgrid(grd)
else:
   sys.exit('wrong format of grid: {}'.format(grd)); sys.stdout.flush()

#load bathymetry on each core
S=zdata(); S.dp=dict(); S.sind=dict()
for m,[fname,psi] in enumerate(zip(fnames,ps)):
    bname=fname.split('.')[0]

    #interpolate depth
    while(True):
        try:
           dpi,sindi=load_dem(gd.x,gd.y,'{}/{}'.format(sdir,fname),fmt=1,position=psi)
           break
        except:
            time.sleep(15)

    #save results
    S.dp[bname]=dpi; S.sind[bname]=sindi
    print('finished reading {},: {}, myrank={}'.format(fname,inum[m],myrank)); sys.stdout.flush()
savez('.load_dem_{}'.format(myrank),S)

#combine results
if ibatch==1: comm.Barrier()
if myrank==0:
   #combine
   S=zdata(); S.dp=dict(); S.sind=dict()
   for i in arange(nproc):
       sname='.load_dem_{}.npz'.format(i)
       Si=read(sname); os.remove(sname)
       S.dp={**S.dp,**Si.dp}
       S.sind={**S.sind,**Si.sind}

   #load bathymetry
   did=zeros(gd.np,'int'); dname=[]
   for i,fname in enumerate(fnames_sort):
       bname=fname.split('.')[0]
       sind=S.sind[bname]; dp=S.dp[bname]
       if reverse_sign==1: dp=-dp #reverse depth sign
       gd.dp[sind]=dp+datum_shift; did[sind]=i+1
       dnamei=[k for k in fnames0 if k.startswith(fname)][0]; dname.append(dnamei) 

   #applying minimum depth
   if regions is not None:
      if len(regions)!=len(rvalues): sys.exit('differet size: regions, rvalues') 
      for i, region in enumerate(regions):
          if not os.path.exists(region): continue
          depth_min=rvalues[i]
          bp=read_schism_bpfile(region,fmt=1)
          sind=inside_polygon(c_[gd.x,gd.y], bp.x,bp.y)
          fp=(sind==1)*(gd.dp<depth_min); gd.dp[fp]=depth_min
          print('finished applying min depth={}: {}'.format(depth_min,region)); sys.stdout.flush()

   #save grid
   if grdout.endswith('npz'): 
      S=zdata(); S.hgrid=gd; savez(grdout,S)
   else:
      gd.write_hgrid(grdout)
   fid=open('{}_dem_id'.format(grdout),'w+'); [fid.write('{}\n'.format(i)) for i in did]; fid.close()
   fid=open('{}_dem_name'.format(grdout),'w+'); [fid.write('{}: {}\n'.format(i+1,k)) for i,k in enumerate(dname)]; fid.close()
   
#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
if ibatch==1: comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora','levante'] else os._exit(0)
