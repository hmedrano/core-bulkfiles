"""
 Script que prepara los datos crudos de FNL y GFS_HD para el uso como 
 forzamientos meteorologicos en la simulacion con NEMO-OPA.

 Se interpolan los datos GFS_HD a la malla de FNL para trabajar con ellos 
 como un solo dataset. 

 El nombre de las variables a procesar se obtiene del archivo de configuracion
 gfsconfig.cfg.

 By Favio Medrano

"""

import os 
import glob
import numpy as np
import netCDF4 as nc 
from scipy import interpolate 
import logging as log
import datetime as dt
# Own libs
import nemoForcingMaker


def findFNL_GFS(searchPath):
    """
     Script que busca en el directorio 'searchPath' por archivos crudos FNL y GFS_HD
     para usarlos como parametro en el metodo doFNL_GFSForcing.
    """
    if os.path.exists(searchPath):
        fnlFile = glob.glob(os.path.join(searchPath,'crudosFNL*.nc'))
        gfsFile = glob.glob(os.path.join(searchPath,'crudosGFS_HD*.nc'))
        if len(fnlFile) > 0 and len(gfsFile) > 0: 
            log.info('FNL File : ' + fnlFile[0])
            log.info('GFS_HD File ' + gfsFile[0])
            return doFNL_GFSForcing(fnlFile[0],gfsFile[0]) 
        else: 
            log.error('No se encontro fnl o gfs en el directorio: ' + searchPath)
            return -1
    else:
        log.error('El directorio ' + searchPath + ' no existe.')
        return -1 


def doFNL_GFSForcing(fnlCrudos,gfsCrudos):
    """
     Metodo que se encarga de leer que variables de los datasets FNL y GFS se van a utilizar como 
     forzamientos.
     Interpola las variables del dataset GFS_HD al la malla de FNL para su manejo como un solo dataset.
     Utiliza el metodo de la clase nemoForcingMaker para generar los archivos de forzamientos.
    """
    
    # Asegurarnos que los archivos fnlCrudos, gfsCrudos y gfsconfig.cfg existan
    if not (os.path.exists(fnlCrudos) and os.path.exists(gfsCrudos) and os.path.exists('gfsconfig.cfg') ):
        log.error('Alguno de los archivos necesarios no existe en el directorio de trabajo: fnlCrudos, gfsCrudos o gfsconfig.cfg')
        return -1
    
    # Leemos el archivo de configuracion, para reconocer que variables 
    # se van a interpolar y preparar para los forzamientos.
    confData = nemoForcingMaker.gfsConfig()
    
    # Leemos los datos e interpolamos GFS_HD a la malla de FNL.
    fnlData = nc.Dataset(fnlCrudos,'r')
    gfsData = nc.Dataset(gfsCrudos,'r') 
    
    
    # interpolar datos de gfsCrudos a la malla de fnl
    # malla gfs 1/2 grado - malla fnl 1/5 grado
    xxn = fnlData.variables['lon'][:]
    yyn = fnlData.variables['lat'][:]
    
    xx = gfsData.variables['lon'][:] 
    yy = gfsData.variables['lat'][:] 
    
    # Variable temporal de gfs y fnl, concatenadas. GFS viene @3hrs FNL @6hrs 
    # Concatenamos @6hrs fnl + gfs_hd
    timeVarFNL = fnlData.variables['time'][:] 
    timeVarGFS = gfsData.variables['time'][:]
    timeFull = np.concatenate((timeVarFNL,timeVarGFS[0:timeVarGFS.size-1:2]))
    
    # Que variables vamos a interpolar:
    lVars = confData.getConfigValueVL('vars')
    newVars = {}
    # Ciclo para interpolar todas las variables
    for var in lVars:
        
        newVars[var] = np.zeros((timeFull.size , yyn.size , xxn.size))  
        if var != 'snodsfc':
            n = 0 
            log.info('Llenando FNL variable: ' + var)
            # Lllenar de datos de FNL 
            for t in range(0,timeVarFNL.size):
                varData = fnlData.variables[var][t][:][:] 
                newVars[var][n][:][:] = varData 
                n = n + 1 
                log.info('N : ' + str(n-1))
                 
            # Llenar datos de GFS, con interpolacion
            log.info('Interpolado GFS variable : ' + var)
            for t in range(0,timeVarGFS.size-1,2):
                log.info('Tiempo ' + str(t))
            
                varData = gfsData.variables[var][t][:][:]
                if np.unique(varData).size == 1:
                    varData = gfsData.variables[var][t+1][:][:]
                    log.info('Tiempo vacio, se tomara el tiempo siguiente!!')
                
                funcInterpol = interpolate.RectBivariateSpline(yy,xx,varData)  
                varDatanew =  funcInterpol(yyn,xxn)
                newVars[var][n][:][:] = varDatanew 
                n = n + 1
                log.info('N : ' + str(n-1))
        else:
            log.info('Dejamos snodsfc con ceros.')
            
    # Cerrar archivos de datos crudos
    fnlData.close()
    gfsData.close()
    
    # Utilizar los scripts para generar archivos mensuales o anuales de los forzamientos
    myForc = nemoForcingMaker.nemoForcing() 
    myForc.makeForcingCoreBulk( {'time' : timeFull, 'lat' : yyn, 'lon': xxn}, newVars, 'yearly' )
    
    return 0    

def selDRange(dst, dFrom, dTo): 
    """ 
     Return num time indices that fit in the "from"(datetime) to "to" parameters in the netCDF4 dataset dst 
    """

    timeVarName = 'time'
    timeV = dst.variables[timeVarName][:]
    fromNum = nc.date2num(dFrom,dst.variables[timeVarName].units) 
    toNum = nc.date2num(dTo, dst.variables[timeVarName].units) 

    return np.argwhere((timeV >= fromNum) & (timeV < toNum)).flatten() 

def getDFile(rawDPath,pDate, dataWildC):
    dtPathW = os.path.join(rawDPath, pDate.strftime('%Y%m%d') , dataWildC) 
    dtFPath = glob.glob(dtPathW) 
    if len(dtFPath) > 0:
        return dtFPath[0]
    else: 
        return False 

def doGFScore_bulk(rawDPath, dataWildC, pivotDate, hdays): 
    # Asegurarnos que los archivos fnlCrudos, gfsCrudos y gfsconfig.cfg existan
    if not (os.path.exists('gfsconfig.cfg') ):
        log.error('Alguno de los archivos necesarios no existe en el directorio de trabajo: fnlCrudos, gfsCrudos o gfsconfig.cfg')
        return -1
    
    # Leemos el archivo de configuracion, para reconocer que variables 
    # se van a interpolar y preparar para los forzamientos.
    confData = nemoForcingMaker.gfsConfig()    
    dtFile = getDFile(rawDPath, pivotDate, dataWildC)
    if (dtFile):
        dst = nc.Dataset(dtFile,'r') 
        yyn = dst.variables['lat'][:]
        xxn = dst.variables['lon'][:] 
        timeSize = dst.variables['time'][:].size
        dst.close()
    else:
        log.error('El archivo con el dataset para la fecha ' + str(pivotDate) + ' no se encontro. Abortando')
        return -1 


    # time variable

    # Que variables procesar.
    lVars = confData.getConfigValueVL('vars') 
    newVars = {}
    # Iniciar el tamano de los arreglos.
    for var in lVars:
        newVars[var] = np.zeros((timeSize + (hdays * 8) , yyn.size , xxn.size ))
    # Time variable
    timeFull = np.zeros((timeSize + (hdays * 8)))
    log.info('Buffer para variables con tamano: ' + str(timeFull.size) )
    log.info('Malla 2D shape: ' + str(yyn.size ) + ' , ' + str(xxn.size )  )


    nI = 0
    # Iniciar los arreglos con los hdays
    for d in range(hdays,0,-1): 
        dtC = pivotDate - dt.timedelta(days=d) 
        dtFPath = getDFile(rawDPath, dtC, dataWildC)
        if (dtFPath):   
            log.info('Obteniendo datos de archivo : ' + str(dtFPath))                 
            dst = nc.Dataset(dtFPath, 'r') 

            # Seleccionar un dia del dataset historico
            dInd = selDRange(dst, dtC , dtC+dt.timedelta(days=1)) 
            for i in dInd:
                for var in lVars:
                    newVars[var][nI][:][:] = dst.variables[var][dInd][:][:]
                timeFull[nI] = dst.variables['time'][dInd]
                nI = nI + 1 

            dst.close()
        else:
            log.info('No se encontro archivo para datos con fecha: ' + str(dtC))


    # Agregar la informacion del pronostico que contenga el dataset "pivotDate" 
    # 
    dtFile = getDFile(rawDPath, pivotDate, dataWildC)
    if (dtFile):
        dst = nc.Dataset(dtFile,'r')
        log.info('Obteniendo datos de archivo : ' + str(dtFile))

        for tI in range(0,timeSize): 
            for var in lVars:
                newVars[var][nI][:][:] = dst.variables[var][tI][:][:] 
            timeFull[nI] = dst.variables['time'][tI] 
            nI = nI + 1 

        dst.close()


    # Utilizar los scripts para generar archivos mensuales o anuales de los forzamientos
    myForc = nemoForcingMaker.nemoForcing() 
    out = myForc.makeForcingCoreBulk( {'time' : timeFull, 'lat' : yyn, 'lon': xxn}, newVars, 3 , 'yearly' )

    return out




def main():
    # Test main.
    #doFNL_GFSForcing('crudosFNL_2014-04-22__2014-04-26.nc','crudosGFS_HD_2014-04-27_00z.nc')
    log.getLogger().setLevel(20)
    # findFNL_GFS('.')
    dd = dt.datetime(2015,1,27)
    doGFScore_bulk('/LUSTRE/hmedrano/STOCK/FORCING-RAW/GFS_RAW','*0P25*.nc', dd , 0)



if __name__ == "__main__":
    main()