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

def main():
    # Test main.
    #doFNL_GFSForcing('crudosFNL_2014-04-22__2014-04-26.nc','crudosGFS_HD_2014-04-27_00z.nc')
    log.getLogger().setLevel(10)
    findFNL_GFS('.')

if __name__ == "__main__":
    main()