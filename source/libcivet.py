import sys
import os
import string
import subprocess
import numpy as np
from scipy import interpolate
import netCDF4 as netcdf4

dir_path = os.path.dirname(os.path.realpath(__file__))
data_root=dir_path+'/../ancillary_data/'

#--  Begin function definitions, labeled as ##(number)
##1
def read_etopo_nc(lon,lat):

    topofile=data_root + 'ETOPO5.nc'

    f1 =  netcdf4.Dataset(topofile, 'r', format='NETCDF4')
    lon1  = np.copy(f1.variables['longitude'][:])
    lat1  = np.copy(f1.variables['latitude'][:])
    elev1 = np.copy(f1.variables['elevation'][:,:])
    f1.close

#    if np.min(lon) < 0:
#        rolldata(elev1,lon1,roll=-1,add360=-1)

    x = interpolate.interp2d(lon1, lat1, elev1, kind='linear')
    elev =x(lon,lat)
  
    minelev = 0.
    ind=[elev < minelev]
    elev[ind]=minelev

    return elev

##2
def plot_coastlines_nc(PlotDataItem,center='Greenwich',pen=(0,0,0),
                       xbounds=None,ybounds=None):

# read in coastline data
    coast_root = data_root+'coasts/'

    coastfile = coast_root + 'coasts_0.20.nc'
    coastfile = coast_root + 'coasts_0.12.nc'

    f1 =  netcdf4.Dataset(coastfile, 'r', format='NETCDF4')
    slength = np.copy(f1.variables['segment_length'][:])
    scoords = np.copy(f1.variables['segment_coordinates'][:,:])
    f1.close

    nm=scoords.shape[0]
    ns=slength.size
    #print 'nm,ns: ',nm,ns

    # for each segment, eliminate those outside bounds and
    # update segment length
    if np.any(xbounds):
        cond=(scoords[:,0] < 0.)
        scoordsp = np.copy(scoords)
        scoordsp[cond,0]+=360
        slengtho=np.copy(slength)
        ind=np.asarray([],dtype=np.int)
        i1 = 0
        for n in range(ns):
            i2=i1+slength[n]
            ind2=np.asarray(range(i1,i2),dtype=np.int)
            i1=np.sum(slength[0:n+1])
            ind3=np.where(np.logical_and(scoordsp[ind2,0] >= xbounds[0],
                                         scoordsp[ind2,0] <= xbounds[1]))[0]
            if ind3.size > 0:
                ind = np.append(ind,ind2[ind3])
            slengtho[n] = ind3.size
                
        scoords=scoords[ind,:]
        ind=np.where(slengtho > 0)[0]
        slength = slengtho[ind]
        nm=scoords.shape[0]
        ns=slength.size
        #print 'nm,ns: ',nm,ns

    if np.any(ybounds):
        slengtho=np.copy(slength)
        ind=np.asarray([],dtype=np.int)
        i1 = 0
        for n in range(ns):
            i2=i1+slength[n]
            ind2=np.asarray(range(i1,i2),dtype=np.int)
            i1=np.sum(slength[0:n+1])
            ind3=np.where(np.logical_and(scoords[ind2,1] >= ybounds[0],
                                         scoords[ind2,1] <= ybounds[1]))[0]
            if ind3.size > 0:
                ind = np.append(ind,ind2[ind3])
            slengtho[n] = ind3.size
                
        scoords=scoords[ind,:]
        ind=np.where(slengtho > 0)[0]
        slength = slengtho[ind]
        nm=scoords.shape[0]
        ns=slength.size
        #print 'nm,ns: ',nm,ns

    # create filter to avoid plotting unconnected segments
    psegs=np.ones(scoords.shape[0],dtype=np.bool_)
    i1=0
    for n in range(ns):
            i2=i1+slength[n]
            psegs[i2-1] = False
            i1=np.sum(slength[0:n+1])
    
    # filter out points that cross greenwich
    if center == 'Dateline':
        for n in range(nm-1):
            if np.logical_and(scoords[n,0] >= 0., scoords[n+1,0] < 0.):
                psegs[n] = False
            if np.logical_and(scoords[n,0] < 0., scoords[n+1,0] >= 0.):
                psegs[n] = False

        cond=(scoords[:,0] < 0.)
        scoords[cond,0]+=360

    # filter out points that cross dateline
    if center == 'Greenwich':
        dthresh = 300.0
        for n in range(nm-1):
            if np.logical_or(np.abs(scoords[n,0]-scoords[n+1,0]) > dthresh,
                             np.abs(scoords[n+1,0]-scoords[n,0]) > dthresh):
                psegs[n] = False

    PlotDataItem.setData(scoords,pen=pen,connect=psegs)

##3
def plot_countries_nc(PlotDataItem,pen=(0,0,0),center='Greenwich',
                      xbounds=None,ybounds=None):
# read in country boundary data
    country_root=data_root+'political_boundaries/'

    countryfile = country_root + 'world_countries.nc'

    f1 =  netcdf4.Dataset(countryfile, 'r', format='NETCDF4')
    slength = np.copy(f1.variables['segment_length'][:])
    scoords = np.copy(f1.variables['segment_coordinates'][:,:])
    f1.close

    nm=scoords.shape[0]
    ns=slength.size
    #print 'nm,ns: ',nm,ns

    # for each segment, eliminate those outside bounds and
    # update segment length
    if np.any(xbounds):
        cond=(scoords[:,0] < 0.)
        scoordsp = np.copy(scoords)
        scoordsp[cond,0]+=360
        slengtho=np.copy(slength)
        ind=np.asarray([],dtype=np.int)
        i1 = 0
        for n in range(ns):
            i2=i1+slength[n]
            ind2=np.asarray(range(i1,i2),dtype=np.int)
            i1=np.sum(slength[0:n+1])
            ind3=np.where(np.logical_and(scoordsp[ind2,0] >= xbounds[0],
                                         scoordsp[ind2,0] <= xbounds[1]))[0]
            if ind3.size > 0:
                ind = np.append(ind,ind2[ind3])
            slengtho[n] = ind3.size
                
        scoords=scoords[ind,:]
        ind=np.where(slengtho > 0)[0]
        slength = slengtho[ind]
        nm=scoords.shape[0]
        ns=slength.size
        #print 'nm,ns: ',nm,ns

    if np.any(ybounds):
        slengtho=np.copy(slength)
        ind=np.asarray([],dtype=np.int)
        i1 = 0
        for n in range(ns):
            i2=i1+slength[n]
            ind2=np.asarray(range(i1,i2),dtype=np.int)
            i1=np.sum(slength[0:n+1])
            ind3=np.where(np.logical_and(scoords[ind2,1] >= ybounds[0],
                                         scoords[ind2,1] <= ybounds[1]))[0]
            if ind3.size > 0:
                ind = np.append(ind,ind2[ind3])
            slengtho[n] = ind3.size
                
        scoords=scoords[ind,:]
        ind=np.where(slengtho > 0)[0]
        slength = slengtho[ind]
        nm=scoords.shape[0]
        ns=slength.size
        #print 'nm,ns: ',nm,ns

    # create filter to avoid plotting unconnected segments
    psegs=np.ones(scoords.shape[0],dtype=np.bool_)
    i1=0
    for n in range(ns):
            i2=i1+slength[n]
            psegs[i2-1] = False
            i1=np.sum(slength[0:n+1])

    # filter out points that cross greenwich
    if center == 'Dateline':
        for n in range(nm-1):
            if np.logical_and(scoords[n,0] >= 0., scoords[n+1,0] < 0.):
                psegs[n] = False
            if np.logical_and(scoords[n,0] < 0., scoords[n+1,0] >= 0.):
                psegs[n] = False

        cond=(scoords[:,0] < 0.)
        scoords[cond,0]+=360

    # filter out points that cross dateline
    if center == 'Greenwich':
        dthresh = 300.0
        for n in range(nm-1):
            if np.logical_or(np.abs(scoords[n,0]-scoords[n+1,0]) > dthresh,
                             np.abs(scoords[n+1,0]-scoords[n,0]) > dthresh):
                psegs[n] = False

    PlotDataItem.setData(scoords,connect=psegs,pen=pen)

##4
def plot_states_nc(PlotDataItem,pen=(0,0,0),center='Greenwich',
                      xbounds=None,ybounds=None):

# read in us state boundary data
    state_root = data_root+'political_boundaries/'

    statefile = state_root + 'us_states.nc'

    f1 =  netcdf4.Dataset(statefile, 'r', format='NETCDF4')
    slength = np.copy(f1.variables['segment_length'][:])
    scoords = np.copy(f1.variables['segment_coordinates'][:,:])
    f1.close

    nm=scoords.shape[0]
    ns=slength.size
    #print 'nm,ns: ',nm,ns

    # for each segment, eliminate those outside bounds and
    # update segment length
    if np.any(xbounds):
        cond=(scoords[:,0] < 0.)
        scoordsp = np.copy(scoords)
        scoordsp[cond,0]+=360
        slengtho=np.copy(slength)
        ind=np.asarray([],dtype=np.int)
        i1 = 0
        for n in range(ns):
            i2=i1+slength[n]
            ind2=np.asarray(range(i1,i2),dtype=np.int)
            i1=np.sum(slength[0:n+1])
            ind3=np.where(np.logical_and(scoordsp[ind2,0] >= xbounds[0],
                                         scoordsp[ind2,0] <= xbounds[1]))[0]
            if ind3.size > 0:
                ind = np.append(ind,ind2[ind3])
            slengtho[n] = ind3.size
                
        scoords=scoords[ind,:]
        ind=np.where(slengtho > 0)[0]
        slength = slengtho[ind]
        nm=scoords.shape[0]
        ns=slength.size
        #print 'nm,ns: ',nm,ns

    if np.any(ybounds):
        slengtho=np.copy(slength)
        ind=np.asarray([],dtype=np.int)
        i1 = 0
        for n in range(ns):
            i2=i1+slength[n]
            ind2=np.asarray(range(i1,i2),dtype=np.int)
            i1=np.sum(slength[0:n+1])
            ind3=np.where(np.logical_and(scoords[ind2,1] >= ybounds[0],
                                         scoords[ind2,1] <= ybounds[1]))[0]
            if ind3.size > 0:
                ind = np.append(ind,ind2[ind3])
            slengtho[n] = ind3.size
                
        scoords=scoords[ind,:]
        ind=np.where(slengtho > 0)[0]
        slength = slengtho[ind]
        nm=scoords.shape[0]
        ns=slength.size
        #print 'nm,ns: ',nm,ns

    # create filter to avoid plotting unconnected segments
    psegs=np.ones(scoords.shape[0],dtype=np.bool_)
    i1=0
    for n in range(ns):
            i2=i1+slength[n]
            psegs[i2-1] = False
            i1=np.sum(slength[0:n+1])

    # filter out points that cross greenwich
    if center == 'Dateline':
        for n in range(nm-1):
            if np.logical_and(scoords[n,0] >= 0., scoords[n+1,0] < 0.):
                psegs[n] = False
            if np.logical_and(scoords[n,0] < 0., scoords[n+1,0] >= 0.):
                psegs[n] = False

        cond=(scoords[:,0] < 0.)
        scoords[cond,0]+=360

    # filter out points that cross dateline
    if center == 'Greenwich':
        dthresh = 300.0
        for n in range(nm-1):
            if np.logical_or(np.abs(scoords[n,0]-scoords[n+1,0]) > dthresh,
                             np.abs(scoords[n+1,0]-scoords[n,0]) > dthresh):
                psegs[n] = False

    PlotDataItem.setData(scoords,connect=psegs,pen=pen)

##5
def plot_rivers_nc(PlotDataItem,pen=(0,0,0),center='Greenwich',
                      xbounds=None,ybounds=None):
# read in river data
    river_root=data_root+'rivers/'
    riverfile = river_root + 'world_rivers.res=4.nc'

    f1 =  netcdf4.Dataset(riverfile, 'r', format='NETCDF4')
    slength = np.copy(f1.variables['segment_length'][:])
    scoords = np.copy(f1.variables['segment_coordinates'][:,:])
    f1.close

    nm=scoords.shape[0]
    ns=slength.size

    # for each segment, eliminate those outside bounds and
    # update segment length
    if np.any(xbounds):
        cond=(scoords[:,0] < 0.)
        scoordsp = np.copy(scoords)
        scoordsp[cond,0]+=360
        slengtho=np.copy(slength)
        ind=np.asarray([],dtype=np.int)
        i1 = 0
        for n in range(ns):
            i2=i1+slength[n]
            ind2=np.asarray(range(i1,i2),dtype=np.int)
            i1=np.sum(slength[0:n+1])
            ind3=np.where(np.logical_and(scoordsp[ind2,0] >= xbounds[0],
                                         scoordsp[ind2,0] <= xbounds[1]))[0]
            if ind3.size > 0:
                ind = np.append(ind,ind2[ind3])
            slengtho[n] = ind3.size
                
        scoords=scoords[ind,:]
        ind=np.where(slengtho > 0)[0]
        slength = slengtho[ind]
        nm=scoords.shape[0]
        ns=slength.size
        #print 'nm,ns: ',nm,ns

    if np.any(ybounds):
        slengtho=np.copy(slength)
        ind=np.asarray([],dtype=np.int)
        i1 = 0
        for n in range(ns):
            i2=i1+slength[n]
            ind2=np.asarray(range(i1,i2),dtype=np.int)
            i1=np.sum(slength[0:n+1])
            ind3=np.where(np.logical_and(scoords[ind2,1] >= ybounds[0],
                                         scoords[ind2,1] <= ybounds[1]))[0]
            if ind3.size > 0:
                ind = np.append(ind,ind2[ind3])
            slengtho[n] = ind3.size
                
        scoords=scoords[ind,:]
        ind=np.where(slengtho > 0)[0]
        slength = slengtho[ind]
        nm=scoords.shape[0]
        ns=slength.size
        #print 'nm,ns: ',nm,ns

    # create filter to avoid plotting unconnected segments
    psegs=np.ones(scoords.shape[0],dtype=np.bool_)
    i1=0
    for n in range(ns):
            i2=i1+slength[n]
            psegs[i2-1] = False
            i1=np.sum(slength[0:n+1])

    # filter out points that cross greenwich
    if center == 'Dateline':
        for n in range(nm-1):
            if np.logical_and(scoords[n,0] >= 0., scoords[n+1,0] < 0.):
                psegs[n] = False
            if np.logical_and(scoords[n,0] < 0., scoords[n+1,0] >= 0.):
                psegs[n] = False

        cond=(scoords[:,0] < 0.)
        scoords[cond,0]+=360

    # filter out points that cross dateline
    if center == 'Greenwich':
        dthresh = 300.0
        for n in range(nm-1):
            if np.logical_or(np.abs(scoords[n,0]-scoords[n+1,0]) > dthresh,
                             np.abs(scoords[n+1,0]-scoords[n,0]) > dthresh):
                psegs[n] = False

    PlotDataItem.setData(scoords,connect=psegs,pen=pen)
                    
