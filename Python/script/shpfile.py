#!/usr/bin/env python3
from pylib import *

def read_shapefile_data(fname):
    import shapefile as shp
    with shp.Reader(fname) as C:
        #----read shapefile----------------
        S=npz_data();
        S.nrec=C.numRecords
        S.type=C.shapeTypeName

        #----read pts------------------------------------------------------------------
        #works for pts and polygon, may not work for other geomerty (need update in this case)
        S.xy=[];
        for i in arange(S.nrec):
            xyi=array(C.shape(i).points);
            parti=array(C.shape(i).parts,dtype='int');
            #insert nan for delimiter
            #to get original index: ind=nonzero(isnan(xyi[:,0]))[0]-arange(len(parti));
            S.xy.append(insert(xyi,parti,nan,axis=0))
        S.xy=squeeze(array(S.xy))

        #---read attributes------------------------------------------------------------
        S.attname=array([C.fields[i][0] for i in arange(1,len(C.fields))]);
        stype=array([type(C.record()[m]) for m in S.attname])
        svalue=array(C.records(),dtype='O');
        S.attvalue=array(zeros(len(S.attname))).astype('O')
        for m in arange(len(S.attname)):
            S.attvalue[m]=svalue[:,m].astype(stype[m])
        S.atttype=stype

        #read prj file if exist---
        bdir=os.path.dirname(os.path.abspath(fname));
        bname=os.path.basename(fname).split('.')[0]
        prjname='{}/{}.prj'.format(bdir,bname)
        if os.path.exists(prjname):
            with open(prjname,'r') as fid:
                S.prj=fid.readline().strip()

    return S

def write_shapefile_data(fname,S,float_len=18,float_decimal=8):
    import shapefile as shp

    #---get nrec-----
    if S.type=='POINT':
        if S.xy.dtype==dtype('O'):
            print('S.xy has a dtype="O" for POINT'); sys.exit()
        else:
            nrec=S.xy.shape[0];
    elif S.type=='POLYLINE' or S.type=='POLYGON':
        if S.xy.dtype==dtype('O'):
            nrec=len(S.xy)
        else:
            nrec=1;

    #---check nrec
    if hasattr(S,'nrec'):
        if nrec!=S.nrec:
            print('nrec inconsistent')
            sys.exit()

    #---write shapefile---------
    with shp.Writer(fname) as W:
        W.autoBalance=1;
        #define attributes
        if hasattr(S,'attname'):
            if S.attvalue.ndim==1:
                stype=[type(S.attvalue[0])]
            elif S.attvalue.ndim==2:
                stype=[type(S.attvalue[m][0]) for m in arange(len(S.attname))]
            for m in arange(len(stype)):
                if stype[m] in [np.int,np.int8,np.int16,np.int32,np.int64]:
                    W.field(S.attname[m],'N')
                elif stype[m] in [np.float,np.float16,np.float32,np.float64]:
                    W.field(S.attname[m],'F',float_len,float_decimal)
                elif stype[m] in [np.str0,np.str,np.str_,np.string_]:
                    W.field(S.attname[m],'C',100)
                else:
                    print('attribute type not included: add here')
                    sys.exit()
        else:
            W.field('field','C')
            W.record()

        #put values
        for i in arange(nrec):
            if S.type=='POINT': #point, W.multipoint(S.xy) is multiple pts features
                vali=S.xy[i]
                W.point(*vali)
            elif S.type=='POLYLINE':
                if S.xy.dtype==dtype('O'):
                    vali=S.xy[i]
                else:
                    vali=S.xy
                #reorganize the shape of vali
                valii=delete_shapefile_nan(vali,0)
                W.line(valii)
            elif S.type=='POLYGON':
                if S.xy.dtype==dtype('O'):
                    vali=S.xy[i]
                else:
                    vali=S.xy
                #reorganize the shape of vali
                valii=delete_shapefile_nan(vali,1)
                W.poly(valii)

            #add attribute
            if hasattr(S,'attname'):
                if S.attvalue.ndim==1:
                    atti=[S.attvalue[i]]
                elif S.attvalue.ndim==2:
                    atti=[S.attvalue[m][i] for m in arange(len(stype))]
                W.record(*atti)

        #----write projection------------
        bname=os.path.basename(fname).split('.')[0]
        bdir=os.path.dirname(os.path.abspath(fname));
        if hasattr(S,'prj'):
            with open('{}/{}.prj'.format(bdir,bname),'w+') as fid:
                fid.write(S.prj)

def delete_shapefile_nan(xi,iloop=0):
    #----delete nan (head and tail), and get ind for the rest
    if xi.ndim==1:
        i1=0; i2=xi.shape[0]
        if isnan(xi[0]): i1=1
        if isnan(xi[-1]): i2=i2-1
        yi=xi[i1:i2]; ind=nonzero(isnan(yi))[0]
    elif xi.ndim==2:
        i1=0; i2=xi.shape[0]
        if isnan(xi[0,0]): i1=1
        if isnan(xi[-1,0]): i2=i2-1
        yi=xi[i1:i2]; ind=nonzero(isnan(yi[:,0]))[0]

    #------reorganize-----------
    if len(ind)==0:
        #close the geomety
        if iloop==1: yi=close_data_loop(yi)

        vi=[yi];
    else:
        vi=[];
        yii=yi[:ind[0]];
        if iloop==1: yii=close_data_loop(yii)
        vi.append(yii)
        for m in arange(len(ind)-1):
            i1=ind[m]+1; i2=ind[m+1];
            yii=yi[i1:i2];
            if iloop==1: yii=close_data_loop(yii);
            vi.append(yii)
        yii=yi[(ind[-1]+1):];
        if iloop==1: yii=close_data_loop(yii)
        vi.append(yii)

    return vi


if __name__=="__main__":
    pass
#    import shapefile as shp
#
#    #---read grid-----
#    gd=read_schism_hgrid('./hgrid.gr3')

#    #---grid bnd--------------------
#    S=npz_data()
#    S.type='POLYLINE'
#    S.nrec=1;
#    for i in arange(gd.nob):
#        ind=gd.iobn[i]
#        xyi=c_[gd.x[ind],gd.y[ind]];
#        xyi=insert(xyi,0,nan,axis=0);
#        if i==0:
#            xy=xyi
#        else:
#            xy=r_[xy,xyi]
#    for i in arange(gd.nlb):
#        ind=gd.ilbn[i]
#        xyi=c_[gd.x[ind],gd.y[ind]];
#        if gd.island[i]==1: xyi=close_data_loop(xyi)
#        xyi=insert(xyi,0,nan,axis=0)
#        xy=r_[xy,xyi]
#    S.xy=xy
#    S.prj=get_prj_file('epsg:26918')

#    #---grid points---------
#    S=npz_data()
#    S.type='POINT'
#    S.xy=c_[gd.x,gd.y]
#    S.prj=get_prj_file('epsg:26918')
#    S.attname=['node_number']
#    S.attvalue=arange(gd.np)+1;
#
#
#
#     #--grid element----------
#    S=npz_data()
#    S.type='POLYGON'
#    elnode=gd.elnode; fp=elnode[:,-1]<0; elnode[fp,-1]=elnode[fp,0]
#    elnode=fliplr(elnode)
#    for i in arange(4):
#        xyi=c_[gd.x[elnode[:,i]],gd.y[elnode[:,i]]]
#        if i==0:
#            xy=xyi[:,:,None]
#        else:
#            xy=c_[xy,xyi[:,:,None]]
#    xy=transpose(xy,[0,2,1]);
#    S.xy=zeros(gd.ne).astype('O')
#    for i in arange(gd.ne):
#        S.xy[i]=xy[i]
#
#    S.attname=['element_number']
#    S.attvalue=arange(gd.ne)+1;
#    S.prj=get_prj_file('epsg:26918')
#
#    write_shapefile_data('test7',S)
#    S0=read_shapefile_data('test7')
