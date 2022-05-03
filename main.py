#
#
# ReadInSARTimeSeries.py : read INSAR time series of cumulation displacement
#                resuling from SBAS Mintpy 
#
# P.Santitamnon (phisan.chula@gmail.com)
# History : v.0.1 , 25 Apr 2022 Initial
#
#
import rasterio as rio 
import pandas as pd
import numpy as np
import h5py
from pathlib import Path
from rasterio.transform import Affine
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import datetime
import json 
from fastapi import FastAPI

from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

origins = [
    "http://127.0.0.1",
    "http://localhost",
    "http://127.0.0.1:5500",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


#############################################################################
class INSAR_TimeSeries:
    def __init__(self, DIR_VELOC ):
        self.DIR_VELOC = DIR_VELOC
        #self.FILE_VELOC   = DIR_VELO.joinpath( 'velocity.h5' )
        #with rio.open( self.FILE_VELOC) as src:
        #    self.VELOC = src.read()
        self.FILE_TIMESR = DIR_VELO.joinpath( 'timeseries.h5' )
        self.ReadH5PY( self.FILE_TIMESR )
        assert( self.ARRAY.shape[0] == len(self.TIMESR) )

    def ReadH5PY( self, FILE_HDF5 ):
        with h5py.File( FILE_HDF5 ) as f:
            #print( f.keys() ); print( f.attrs.keys() )
            self.LENGTH     = f.attrs['LENGTH']
            self.WIDTH      = f.attrs['WIDTH']
            self.EPSG       = f.attrs['EPSG']
            self.UTM_ZONE   = f.attrs['UTM_ZONE']
            self.ORBIT_DIR  = f.attrs['ORBIT_DIRECTION']
            self.START_DATE = f.attrs['START_DATE']
            self.END_DATE   = f.attrs['END_DATE']
            self.NUM_IFG    = int( f.attrs['mintpy.networkInversion.numIfgram'] )
            ################################# 
            self.X_FIRST = float( f.attrs['X_FIRST'] )
            self.X_STEP  = float( f.attrs['X_STEP']  )
            self.Y_FIRST = float( f.attrs['Y_FIRST'] )
            self.Y_STEP  = float( f.attrs['Y_STEP']  )
            #################################
            self.TR = Affine( self.X_STEP, 0.0         , self.X_FIRST, 
                                      0.0, self.Y_STEP , self.Y_FIRST  )
            timesr = pd.DataFrame( list(f['date']), columns=['Date'] ) 
            timesr['Date'] = pd.to_datetime( timesr['Date'].str.decode('utf-8') )
            timesr['DayCnt'] = timesr['Date']-timesr.iloc[0]['Date']
            timesr['DayCnt'] = timesr.DayCnt.dt.days
            self.TIMESR = timesr
            self.ARRAY = np.stack( f['timeseries'] )
   
    def Tr_XY(self, row_col ):
        return self.TR*row_col

    def Tr_RowCol(self, xy ):
        return ~self.TR*xy

    def getTimeSeriesPnt( self, xy ):
        r,c = self.Tr_RowCol( xy ) 
        df = self.TIMESR.copy()
        df['CumDispl_m'] = self.ARRAY[ : , int(c), int(r) ]  # cm
        coeff = np.polyfit( df.DayCnt, df.CumDispl_m, 1 )  # linear regress
        #import pdb; pdb.set_trace()
        def EstDisp( row, coeff ):
            return row.DayCnt*coeff[0] + coeff[1]
        df['CumDispl_m_'] = df.apply( EstDisp, axis='columns', args=[coeff] ) 
        velo_cmy = 100*coeff[0]*365  # average cm per year
        return df,velo_cmy
            
    def PlotTimeSeries( self, dfTIMESR, TITLE, PLOT_FILE=None ):
        fig,ax = plt.subplots()
        dfTIMESR.plot.scatter( "Date", "CumDispl_m", alpha=0.5, 
                               color='r', rot=45, ax=ax )
        dfTIMESR.plot.line( "Date", "CumDispl_m_", alpha=0.5, 
                            lw=5,   color='b', rot=45, ax=ax )
        ax.set_xlabel( '' )
        ax.set_ylabel( 'Cumulation Displacement (m)')
        ax.set_title( TITLE )
        plt.tight_layout()
        if PLOT_FILE is None:
            plt.show()
        else:
            plt.savefig( PLOT_FILE )


#####################################################
#Band 1: -0.0340674
pnt = { 'PLUTALUANG' :  ( 711_666, 1_405_782 ),  }   #  -0.0340674  mm/year
DIR_VELO = Path('')

insar = INSAR_TimeSeries( DIR_VELO )
list_x=[]
list_y=[]
for k,v in pnt.items():
    ts,velo_cmy = insar.getTimeSeriesPnt( v  )
    list_ts = ts.values.tolist()
    for i in list_ts:
      list_x.append(i[0].strftime("%Y-%m-%d"))
      list_y.append(i[2])
    # title = f'{DIR_VELO} "{k}" Avg.Velo.:{velo_cmy:+.2f} cm/year'
    # PLOT_FILE = f'{k}.png' 
    # print( f'Plotting {PLOT_FILE}...' )
    # insar.PlotTimeSeries( ts , title, PLOT_FILE )
json_result = {'CumDispl_m_':list_y,'date':list_x}

@app.get("/")
def root():
    return json.dumps(json_result)

#import pdb; pdb.set_trace()
