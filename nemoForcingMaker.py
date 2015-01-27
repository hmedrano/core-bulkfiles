"""
 Scripts para generar archivos de forzamientos meteorologicos para el proceso de 
 simulacion con NEMO-OPA 

 By Favio Medrano 

  2014-08-11 : El metodo makeForcingCoreBulk produce los archivos de forzamientos para 
               la simulacion con el NEMO-OPA, el formato de salida son en archivos mensuales
               o anuales.
             : La variable de dimension temporal que se le agrega a los archivos de forzamientos
               se genera con el metodo dateToNemoCalendar. 

"""

import os
import logging as log
from ConfigParser import ConfigParser
import netCDF4 as nc 
import numpy as np 
from scipy import interpolate 
import datetime as dt
# own libs
import netcdfFile


class gfsConfig:
        """
         Clase padre, que se encarga de cargar el archivo de configuracion gfsconfig.cfg
         Contiene metodos para acceder a las entradas del archivo de configuracion.
        """
        configfile = 'gfsconfig.cfg'
        configData = None        

        def __init__(self):
                # Leer archivo de configuracion
                self.readConfig()

        def readConfig(self):
                """
                 Funcion que lee el archivo de configuracion, llamada desde el constructor de la clase.
                """
                self.configData = ConfigParser()
                if os.path.exists(self.configfile):
                        self.configData.read(self.configfile)
                        return 1
                else:
                        log.warning('Archivo de configuracion no existe! (' + self.configfile + ')')
                        return -1
                        
        def getKeyValue(self,key,value):
                try:
                    if self.configData != None:
                        rvalue = self.configData.get(key,value)
                        log.info('Llave : ' + value + ' Atributo leido: ' + str(rvalue))
                        return rvalue
                    else:
                        log.warning('Configuracion: ' + value + ', no existe en el archivo de configuracion!')
                        return ''
                except:
                    return ''
                                             
        def getConfigValue(self,value):
                """ 
                Devuelve la los datos de la llave "value", para el grupo 'gfs_data'    
                """  
                return self.getKeyValue('gfs_data',value)
             
        def getConfigValueV(self,value):
                """ 
                Devuelve la los datos de la llave "value", para el grupo 'variables'    
                """   
                return self.getKeyValue('variables',value)   
        

        def getConfigValueVL(self,value):
                """ 
                Devuelve los datos de la llave "value", que vienen en formato de lista separada por comas
                regresa, un <python list> con los datos.    
                """            
                raw = self.getKeyValue('variables',value)
                return [ e.strip() for e in raw.replace('\n','').split(',') ]
             
        def datasetExists(self,fname):
                try:
                    dst = nc.Dataset(fname,'r') 
                    dst.close() 
                    return True
                except:
                    return False   
             


class nemoForcing(gfsConfig):
    """
     Clase encargada de crear archivo de forzamientos, adaptados para el
     modelo nemo-opa
    """   
    
    variablesRename = {'ugrd10m' : 'u10', 'vgrd10m' : 'v10', 'tcdcclm' : 'tcc' , 'tmp2m' : 't2', 'spfh2m' : 'q2', 'dlwrfsfc' : 'radlw', 'dswrfsfc' : 'radsw' , 'pratesfc' : 'precip', 'snodsfc' : 'snow'}
    
    def datetime2matlabdn(self,fecha):
        mdn = fecha # + timedelta(days = 366)
        frac = (fecha-dt.datetime(fecha.year,fecha.month,fecha.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
        return mdn.toordinal() + frac


    def dateToNemoCalendar(self,data, ctype='gregorian',give='full'):
        """ 
         Codigo tomado del codigo de nemo IOIPSL/src/calendar.f90 para construir el valor de la variable temporal en modo ordinal.
         segun el calendario con que se prepare la configuracion de nemo. Estos pueden ser:
          gregorian, noleap, all_leap, 360_day, julian

         El parametro 'give' se utiliza para escojer el valor que regresa la funcion:
          full (default) : Regresa el valor ordinal al que corresponde la fecha 'data' (datetime) segun
                           el calendario quese haya seleccionado en 'ctype' 
          monthLen       : Regresa los dias que tiene el mes contenido en la fecha 'data' (datetime), segun
                           el calendario que se haya seleccionado en 'ctype'
          yearLen        : Regresa los dias que contiene el ano, segun el calendario que se haya seleccionado
                           en 'ctype'

         Se utiliza como fecha epoch 1-1-1950 
        """
        epochy = 1950
        epochm = 1
        epochd = 1
        if ctype == 'gregorian':
            oneyear = 365.2425
            ml = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        elif ctype == 'noleap':
            oneyear = 365 
            ml = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        elif ctype == 'all_leap': # 366_day
            oneyear = 366.0
            ml = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
        elif ctype == '360_day':
            oneyear = 360.0 
            ml = np.array([30,30,30,30,30,30,30,30,30,30,30,30])
        elif ctype == 'julian':
            oneyear = 365.25
            ml = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        
        if give == 'yearLen':
            return oneyear

        if give == 'monthLen':
            return ml[data.month -1] 

        if (not isinstance(data,np.ndarray)):
            data = np.array([data])

        newc = np.zeros((len(data)),float)
        for idx,v in enumerate(data):
            y = v.year 
            m = v.month - 1
            d = (v.day-1) + (v.hour / 24.0) 
            nnumdate = (y - epochy) * oneyear 
            for nm in range(0,m):
                nnumdate = nnumdate + (ml[nm])
            nnumdate = nnumdate + d 
            newc[idx] = nnumdate 

        return np.squeeze(newc) 

    def mtDToDatetime(self, dataOrd):
        return dt.datetime.fromordinal(int(dataOrd)) + dt.timedelta(days=dataOrd%1) - dt.timedelta(days=1)

    def makeForcingCoreBulk(self,dimsData,varsData,timeD , sFileSize='yearly'):
        """
         dimsData es un <python dict> con el siguiente formato:
          {'time' : values , 'lat' : values , 'lon' : values}
           
         varsData es un <python dict> con el siguiente formato:
          {'var1' : values , 'var2' : values, ...}

         timeD Indica a cada cuantas horas vienen los datos en varsData, puede ser @ 3hrs, 6hrs, 12hrs

         La funcion hace un ciclo con los valores de la dimension temporal, se obtiene el mes y el ano al que pertenece ese instante
         cada vez que el mes, en el siguiente instante de tiempo, se crea un nuevo archivo con datos "mensuales" o "anuales" 
         Se hace una excepcion para la variable radsw, pues estos datos se guardan con una periodicidad diaria, a diferencia de las 
         demas variables que son cada 6hrs. 
        """
        currentMFile = None
        ncFiles = {}
        lVars = self.getConfigValueVL('vars')
        vUnit = self.getConfigValueVL('units')
        vLN = self.getConfigValueVL('longnames')

        sCalendarType = 'noleap'
        
        N = 0 
        RADSWC = 0 
        dfD = ((24 / timeD) - 1)
        # Ciclo en la variable de dimension de tiempo: 
        for idx_tval,tval in enumerate(dimsData['time'][:]):
            # Convertir el valor tval a datetime.
            # El formato ordinal de python es prolectic gregorian, por lo que se le resta un dia para ajustarlo
            # del formato standard gregorian del que vienen los datos en nomads.
            #tvaldate = dt.datetime.fromordinal(int(tval)) + dt.timedelta(days=tval%1) - dt.timedelta(days=1)
            tvaldate = self.mtDToDatetime(tval) 
            tvaldate_year = tvaldate.year 
            tvaldate_month = tvaldate.month
            # 
            # Crear archivos de salida, segun corresponda a la fecha actual.
            # Se crearan archivos mensuales o anuales de datos.
            #
            compareVarDummyForFileSize = tvaldate.year if (sFileSize == 'yearly') else tvaldate.month

            if currentMFile == None or currentMFile != compareVarDummyForFileSize:

                # Crear el archivo mensual o anual segun coresponda
                if sFileSize == 'yearly':
                    currentMFile = tvaldate.year 
                else:
                    currentMFile = tvaldate.month

                # Crear valores de dimension para el actual mes
                # ds es el primer dia del (mes o ano) del archivo currentMFile
                if sFileSize == 'yearly':
                    ds = dt.datetime(tvaldate.year, 1,1) 
                else: 
                    ds = dt.datetime(tvaldate_year,tvaldate_month,1)

                nm = (tvaldate_month + 1) if tvaldate_month < 12 else 1
                ny = tvaldate_year if tvaldate_month < 12 else (tvaldate_year + 1)
                de = dt.datetime(ny,nm,1) 

                # Tamano de la dimension temporal de este periodo, segun el calendario 
                # que se utilize 
                if sFileSize == 'yearly':
                    sFileTDimSize = self.dateToNemoCalendar(tvaldate,sCalendarType,'yearLen')
                else:
                    sFileTDimSize = self.dateToNemoCalendar(tvaldate,sCalendarType,'monthLen')

                # timeVD contiene los valores temporales para el (mes o ano) que se esta trabajando.
                # espaciado cada "timeD" horas
                timeVD=[]
                timeVDD= []
                for t in range(0,(sFileTDimSize*24) / timeD):
                    #timeVD.append( self.datetime2matlabdn( ds + dt.timedelta(hours=6.0*t) + dt.timedelta(days=1) ) ) 
                    timeVD.append( self.dateToNemoCalendar(ds + dt.timedelta(hours=timeD * t),sCalendarType) )

                # timeVDD contiene los valores temporales para el (mes o ano), espaciado en dias 
                for t in range(0,sFileTDimSize):
                    #timeVDD.append( self.datetime2matlabdn( ds + dt.timedelta(days=1*t) + dt.timedelta(days=1) ) )
                    timeVDD.append( self.dateToNemoCalendar(ds + dt.timedelta(days=1.0*t) , sCalendarType) )
                #
                # Crear archivos netcdf a los cuales se descargaran los datos.
                for var in self.variablesRename.keys():
                    log.info('Procesando variable : ' + var)
                    ncFiles[var] = netcdfFile.netcdfFile()
                    if sFileSize == 'yearly':
                        ncFiles[var].createFile('drowned_' + self.variablesRename[var] + '_GFS_y' + str(tvaldate.year) + '.nc' )
                    else:
                        ncFiles[var].createFile('drowned_' + self.variablesRename[var] + '_GFS_y' + str(tvaldate.year) + '_M' + ("%02d"%tvaldate.month) + '.nc' )

                    # Crear dimensiones y sus variables para referencia.
                    ncFiles[var].createDims({'time':None , 'lat' : dimsData['lat'].size , 'lon' : dimsData['lon'].size})
                    dimVars = { 'time' : { 'dimensions': ['time']  , 'attributes' : {'units':'days since 1950-01-01 00:00:00', 'time_origin' : '1950-01-01 00:00:00', 'calendar' : 'noleap'} , 'dataType' : 'f8' }  
                               ,'lat' :  { 'dimensions': ['lat']   , 'attributes' : {'units':'degree_north'} , 'dataType' : 'f8' }  
                               ,'lon' :  { 'dimensions': ['lon']   , 'attributes' : {'units':'degree_east'}  , 'dataType' : 'f8' }  }
                    ncFiles[var].createVars(dimVars)
                    # Salvar datos de variables de dimension.
                    # 
                    if var == 'dswrfsfc':
                        ncFiles[var].saveData({'time':np.array(timeVDD),'lat' : dimsData['lat'], 'lon' : dimsData['lon']})
                    else:
                        ncFiles[var].saveData({'time':np.array(timeVD),'lat' : dimsData['lat'], 'lon' : dimsData['lon']})
                    # Indice de el arreglo units y long names
                    dindex = lVars.index(var)
                    # Crear la definicion de la variable en el archivo.
                    ncFiles[var].createVars({self.variablesRename[var] : {'dimensions' : ['time','lat','lon'] , 'attributes' : {'units' : vUnit[dindex], 'long_name' : vLN[dindex] , '_FillValue' : 9.999e+20 } , 'dataType' : 'f4' } })
             
            # Salvar Dato de variable en su archivo correspondiente.
            for var in self.variablesRename.keys():
                if var != 'dswrfsfc':
                    # Localizar en que indice salvar el dato: 
                    idx = (np.abs(timeVD - self.dateToNemoCalendar(tvaldate,sCalendarType))).argmin()
                    log.info('Salvando en el indice ' + str(idx) + ' timeVD: ' + str(timeVD[idx]) + '  tvaldate : ' + str(tvaldate) + '  N=' + str(N+1))
                    #ncFiles[var].saveData({self.variablesRename[var] : varsData})
                    ncFiles[var].saveDataS(self.variablesRename[var],varsData[var][N][:][:],(idx))
                    if idx_tval == 0:
                        for di in range(0,idx):
                            ncFiles[var].saveDataS(self.variablesRename[var],varsData[var][N][:][:],(di))
                    if idx_tval >= (dimsData['time'][:].size - 1):
                        for di in range(idx,len(timeVD)):
                            ncFiles[var].saveDataS(self.variablesRename[var],varsData[var][N][:][:],(di))
                else:
                    if RADSWC >= dfD: 
                        #TODO: Verificar que se estan haciendo promedios diarios, las fechas: dimsData['time'][N-3] .. dimsData['time'][N-3] deben pertenecer
                        #      al mismo dia.
                        idxD = (np.abs(timeVDD - self.dateToNemoCalendar(self.mtDToDatetime(dimsData['time'][N - dfD ]), sCalendarType))).argmin() 
                        dataDaily = 0 
                        fechastr = ''
                        for c in range(0,dfD+1):
                            dataDaily = dataDaily  +  varsData[var][N - c][:][:]
                            fechastr = str(self.dateToNemoCalendar(self.mtDToDatetime(dimsData['time'][ N -c ]),sCalendarType)) + fechastr + ' '
                        dataDaily = dataDaily / (dfD+1)
                        #dataDaily = ( varsData[var][N][:][:] + varsData[var][N-1][:][:] + varsData[var][N-2][:][:] + varsData[var][N-3][:][:] ) / 4.0
                        ncFiles[var].saveDataS(self.variablesRename[var], dataDaily ,(idxD)) 
                        #fechastr = str(self.dateToNemoCalendar(self.mtDToDatetime(dimsData['time'][N-3]),sCalendarType)) + ' ' + \
                        #           str(self.dateToNemoCalendar(self.mtDToDatetime(dimsData['time'][N-2]),sCalendarType)) + ' ' + \
                        #           str(self.dateToNemoCalendar(self.mtDToDatetime(dimsData['time'][N-1]),sCalendarType)) + ' ' + \
                        #           str(self.dateToNemoCalendar(self.mtDToDatetime(dimsData['time'][N]),sCalendarType))  
                        log.info('Haciendo promedio diario para radsw con fechas : ' + fechastr)
                        log.info('Salvando en el indice : ' + str(idxD))

                        # Rellenar registros de datos al inicio y al final del archivo 
                        if idx_tval == dfD:
                            log.info('Rellenando inicio de archivo, radsw')
                            for di in range(0,idxD):
                                ncFiles[var].saveDataS(self.variablesRename[var],dataDaily,(di))
                        if idx_tval >= int(dimsData['time'][:].size / (dfD+1)) * (dfD+1) - 1:
                            log.info('Rellenando final de archivo, radsw')
                            for di in range(idxD,len(timeVDD)):
                                ncFiles[var].saveDataS(self.variablesRename[var],dataDaily,(di))                        
                        
            # Contador indice en variable temporal para datos @"timeD"hrs            
            N = N + 1 
            # Contador indice en variable temporal para "radsw", @24 hrs.
            RADSWC = RADSWC + 1
            if RADSWC > dfD:
                RADSWC = 0
        
        log.info('makeForcingCoreBulk: Informacion salvada.')
        return 0
        
   
   
    
    