# This module provides more interactive plotting environment meant for the jupyter notebook.


import alphatims.bruker
import alphatims.utils
import holoviews as hv
from holoviews import opts
import numpy as np
import datashader as ds
from holoviews.operation.datashader import rasterize, dynspread
from holoviews.operation import decimate
from holoviews import streams
import sqlite3
import pandas as pd


# This class is for plotting 4D data
# abstract class for plotting
class Plotter4D:
    def __init__(self, fIn,
            cmbOpts = dict(tools=['hover']), 
            mobOpts = dict(width=500, height=200, color='blue', ylim=(0, 50000)), 
            heatOpts = dict(width=500, height=500, cmap='inferno',  cnorm='eq_hist', bgcolor='black', colorbar=True, xlim=(0.6, 1.5), ylim=(200, 1800)), 
            specOpts = dict(width=200, height=500, color='blue', invert_axes=True, ylim=(0, 100000)), exp=None, decimateMob=True):

        if exp is None: # load experiment object (can be preloaded to save time
            alphatims.utils.set_threads(4)
            self.exp = alphatims.bruker.TimsTOF(fIn)
        self.conn = sqlite3.connect(fIn + "/analysis.tdf")


        self.exp = exp
        self.cmbOpts = cmbOpts
        self.heatOpts = heatOpts
        self.specOpts = specOpts
        self.mobOpts = mobOpts

        self.frames = [0] # rt frames to plot, just the first one by default

        self.decimateMob = decimateMob

    # Define plotting functions
    def plotHeatMapAcrossRT(self, frame):
        # fetch data
        data = hv.Dataset(self.exp[frame][['mz_values', 'mobility_values', 'intensity_values']], kdims=['mobility_values', 'mz_values'], vdims=['intensity_values'])

        heat = data.to(hv.Points).opts(**self.cmbOpts).opts(**self.heatOpts)
        return heat
    
    def plotSpectrumAcrossRT(self, frame, x_range):
        data = hv.Dataset(self.exp[frame][['mz_values', 'mobility_values', 'intensity_values']], kdims=['mobility_values', 'mz_values'], vdims=['intensity_values'])
        dataSpec = data.redim(intensity_values=hv.Dimension('intensity_across_im'))
        spectrum = dataSpec.aggregate('mz_values', np.sum).to(hv.Spikes).opts(**self.cmbOpts).opts(**self.specOpts)
        return spectrum

    def plotMobilogramAcrossRT(self, frame, y_range):
        data = hv.Dataset(self.exp[frame][['mz_values', 'mobility_values', 'intensity_values']], kdims=['mobility_values', 'mz_values'], vdims=['intensity_values'])
        dataMob = data.redim(intensity_values=hv.Dimension('intensity_across_mz'))
        mobilogram = dataMob.select(mz_values=y_range).aggregate('mobility_values', np.sum).to(hv.Curve).opts(**self.cmbOpts).opts(**self.mobOpts)
        return mobilogram

    def plot(self):
        # heat
        heatRt = hv.DynamicMap(self.plotHeatMapAcrossRT, kdims=['frame']).redim.values(frame=self.frames)
        heatRt = dynspread(rasterize(heatRt, aggregator=ds.sum('intensity_values')), threshold=0.9, shape='square').opts(**self.cmbOpts).opts(**self.heatOpts)

        # spectrum
        rangeMob = streams.RangeX((0.6, 1.5), source=heatRt)
        specRt = hv.DynamicMap(self.plotSpectrumAcrossRT, streams=[rangeMob], kdims=['frame']).redim.values(frame=self.frames).opts(**self.cmbOpts).opts(**self.specOpts)

        #mobilogram
        rangeMz = streams.RangeY((400, 1200), source=heatRt)
        if self.decimateMob:
            print("Decimating ...")
            mobRt = decimate(hv.DynamicMap(self.plotMobilogramAcrossRT, streams=[rangeMz], kdims=['frame']).redim.values(frame=self.frames).opts(**self.cmbOpts).opts(**self.mobOpts))
        else:
            mobRt = hv.DynamicMap(self.plotMobilogramAcrossRT, streams=[rangeMz], kdims=['frame']).redim.values(frame=self.frames).opts(**self.cmbOpts).opts(**self.mobOpts)

        return (heatRt + specRt + mobRt).cols(2)





# this class is for getting an MS1 plotter object
class MS1Plotter(Plotter4D):
    def __init__(self, fIn, **kwargs):
        super().__init__(fIn, **kwargs)

        # add appropriate ranges for MS1 view
        self.mobOpts['ylim'] = (0, 50000)

        self.heatOpts['xlim'] = (0.6, 1.5)
        self.heatOpts['ylim'] = (200, 1800)

        self.specOpts['ylim'] = (0, 100000)
    
    # return list of MS1 frames
    def fetchMS1Frames(self):
        return list(self.exp.frames[self.exp.frames['MsMsType'] == 0]['Id'].values)


    # this method returns the 4D plot
    def plot(self):

        # fetch MS1 frames list
        self.frames = self.fetchMS1Frames()

        return super().plot()


# This plots a swath
class SwathPlotter(Plotter4D):

    def __init__(self, fIn, **kwargs):
        super().__init__(fIn, **kwargs)
    
    def fetchSwathFrames(self, swathGroup):
        query = "select frame from DiaFrameMsMsInfo where windowGroup={}".format(swathGroup)
        rslt = self.conn.execute(query)
        return [i[0] for i in rslt.fetchall()]


    def getSwathMetaInfo(self, swathGroup, multiplex):
        query = "select ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth from DiaFrameMsMsWindows where WindowGroup = {}".format(swathGroup)

        out = pd.read_sql(query, self.conn)



        out['ImBegin'] = self.exp.convert_from_indices(scan_indices = out['ScanNumBegin'].values, return_mobility_values=True)['mobility_values']
        out['ImEnd'] = self.exp.convert_from_indices(scan_indices = out['ScanNumEnd'].values, return_mobility_values=True)['mobility_values']

        return out.iloc[4 - multiplex]

    # swathNum is the group number (usually between 1-16)
    # multiplex should be value between 1 and 4)
    def plot(self, swathNum, multiplex):
        # fetch frame list

        self.frames = self.fetchSwathFrames(swathNum)

        meta = self.getSwathMetaInfo(swathNum, multiplex)


        print("Meta Info On Swath:\n{}".format(meta))


        # set the bounds appropriate for the swath
        self.mobOpts['ylim'] = (0, 20000)

        self.heatOpts['xlim'] = (meta['ImEnd'], meta['ImBegin'])
        self.specOpts['ylim'] = (0, 10000)


        return super().plot()
