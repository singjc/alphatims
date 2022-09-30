# This module provides more interactive plotting environment meant for the jupyter notebook.


import alphatims.bruker
import alphatims.utils
import alphatims.osw_plotter
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
# abstract class for plotting shold not be used
class Plotter4D:
    def __init__(self, fIn,
            cmbOpts = dict(tools=['hover']), 
            mobOpts = dict(width=500, height=200, color='blue', ylim=(0, 50000)), 
            heatOpts = dict(width=500, height=500, cmap='inferno',  cnorm='eq_hist', bgcolor='black', colorbar=True, xlim=(0.6, 1.5), ylim=(200, 1800)), 
            specOpts = dict(width=200, height=500, color='blue', invert_axes=True, ylim=(0, 100000)), 
            chromOpts = dict(width=500, height=200, color='blue', ylim=(0, 10000)),
            pointerChromOpts = dict(width=500, height=200, color='green', ylim=(0, 10000)), exp=None, decimateMob=True):

        if exp is None: # load experiment object (can be preloaded to save time
            alphatims.utils.set_threads(4)
            self.exp = alphatims.bruker.TimsTOF(fIn)
        self.conn = sqlite3.connect(fIn + "/analysis.tdf")


        self.exp = exp
        self.cmbOpts = cmbOpts
        self.heatOpts = heatOpts
        self.specOpts = specOpts
        self.mobOpts = mobOpts
        self.chromOpts = chromOpts
        self.pointerChromOpts = pointerChromOpts


        ### debug variables
        self.mobRange = None
        self.curFrame = None
        self.mzRange = None



        self.frames = [0] # rt frames to plot, just the first one by default

        self.decimateMob = decimateMob

        # update the chromgram options based on supplied frames
        self.chromOpts['xlim'] = (min(self.frames)-1, max(self.frames) + 1)
        self.pointerChromOpts['xlim'] = self.chromOpts['xlim']

    # this function corrcts from where mouse was clicked to the closest valid frame value
    def correctFrame(self, frame):
        absolute_difference_function = lambda list_value : abs(list_value - frame)
        return min(self.frames, key=absolute_difference_function)


    # Define plotting functions
    def plotHeatMapAcrossRT(self, frame, y=0):


        # "correct" frame position to an actual query frame
        closest_frame = self.correctFrame(frame)

        self.curFrame = closest_frame

        # fetch data
        data = hv.Dataset(self.exp[closest_frame][['mz_values', 'mobility_values', 'intensity_values']], kdims=['mobility_values', 'mz_values'], vdims=['intensity_values'])

        heat = data.to(hv.Points).opts(**self.cmbOpts).opts(**self.heatOpts)
        return heat
    
    def plotSpectrumAcrossRT(self, frame, x_range, y=0):
        # "correct" frame position to an actual query frame
        closest_frame = self.correctFrame(frame)

        self.curFrame = closest_frame
        self.mobRange = x_range


        # fetch data
        data = hv.Dataset(self.exp[closest_frame][['mz_values', 'mobility_values', 'intensity_values']], kdims=['mobility_values', 'mz_values'], vdims=['intensity_values'])
        dataSpec = data.redim(intensity_values=hv.Dimension('intensity_across_im'))

        # only select the spectrum that are currently visible in theat heat map and then aggregate
        spectrum = dataSpec.select(mobility_values=x_range).aggregate('mz_values', np.sum).to(hv.Spikes).opts(**self.cmbOpts).opts(**self.specOpts)
        return spectrum

    def plotMobilogramAcrossRT(self, frame, y_range, y=0):
        # "correct" frame position by mouse to an actual query frame
        closest_frame = self.correctFrame(frame)

        self.curFrame = closest_frame
        self.mzRange = y_range

        data = hv.Dataset(self.exp[closest_frame][['mz_values', 'mobility_values', 'intensity_values']], kdims=['mobility_values', 'mz_values'], vdims=['intensity_values'])
        dataMob = data.redim(intensity_values=hv.Dimension('intensity_across_mz'))
        mobilogram = dataMob.select(mz_values=y_range).aggregate('mobility_values', np.sum).to(hv.Curve).opts(**self.cmbOpts).opts(**self.mobOpts)
        return mobilogram
    
    # for the given range "shown" in the heatmap plot the total intensity
    def plotChromatogram(self, x_range, y_range, y=0):
        data = hv.Dataset(self.exp[self.frames][['mz_values', 'mobility_values', 'intensity_values', 'frame_indices']], kdims=['mobility_values', 'mz_values', 'frame_indices'], vdims=['intensity_values'])
        dataChrom = data.redim(intensity_values=hv.Dimension('total_intensity'))
        chromatogram  = dataChrom.select(mz_values=y_range, mobility_values=x_range).aggregate(['frame_indices'], np.sum).to(hv.Curve).opts(**self.cmbOpts).opts(**self.chromOpts)
        return chromatogram


    ### helper functions for plotChromatogram ###
    ## plot a line corresponding to the selected chromatogram
    def plotChromatogramPointer(self, frame,y=0):
        # get the closet point to the x value to the clicked region
        closest_value = self.correctFrame(frame)

        # return graph of vertical line
        return hv.VLine(closest_value).opts(**self.pointerChromOpts).opts(title=f'Frame: {closest_value}')


    def plot(self, plotChrom=False):
        # heat
        heatRt = hv.DynamicMap(self.plotHeatMapAcrossRT, kdims=['frame']).redim.values(frame=self.frames)
        heatRt = dynspread(rasterize(heatRt, aggregator=ds.sum('intensity_values')), threshold=0.9, shape='square').opts(**self.cmbOpts).opts(**self.heatOpts)

        # spectrum
        rangeMob = streams.RangeX(x_range=(0.6, 1.5), source=heatRt, transient=True)
        specRt = hv.DynamicMap(self.plotSpectrumAcrossRT, streams=[rangeMob], kdims=['frame']).redim.values(frame=self.frames).opts(**self.cmbOpts).opts(**self.specOpts)

        #mobilogram
        rangeMz = streams.RangeY(y_range=(400, 1200), source=heatRt, transient=True)
        if self.decimateMob:
            mobRt = decimate(hv.DynamicMap(self.plotMobilogramAcrossRT, streams=[rangeMz], kdims=['frame']).redim.values(frame=self.frames).opts(**self.cmbOpts).opts(**self.mobOpts))
        else:
            mobRt = hv.DynamicMap(self.plotChromatogram, streams=[rangeMz], kdims=['frame']).redim.values(frame=self.frames).opts(**self.cmbOpts).opts(**self.mobOpts)

        # Chromatogram if selected
        if plotChrom:
            #chrom = hv.DynamicMap(self.plotMobilogramAcrossRT, streams=[rangeMob, rangeMz], kdims=['frame']).redim.values(frame=self.frames).opts(**self.cmbOpts).opts(**self.chromOpts)
            chrom = hv.DynamicMap(self.plotChromatogram, streams=[rangeMz, rangeMob]).opts(**self.cmbOpts).opts(**self.chromOpts)

        return (heatRt + specRt + mobRt + chrom ).cols(2)


    def plotChrom(self):
        # This is a stream of the current frame, can be created before because it is transient, and needed for initial plotting of other plots
        # create 
        framePointer = streams.SingleTap(x=self.frames[0], y=None, rename={'x':'frame'}, transient=True) 

        # heat
        heatRt = hv.DynamicMap(self.plotHeatMapAcrossRT, streams=[framePointer])
        heatRt = dynspread(rasterize(heatRt, aggregator=ds.sum('intensity_values')), threshold=0.9, shape='square').opts(**self.cmbOpts).opts(**self.heatOpts)

        # spectrum
        rangeMob = streams.RangeX(x_range=(0.6, 1.5), source=heatRt)
        specRt = hv.DynamicMap(self.plotSpectrumAcrossRT, streams=[rangeMob, framePointer]).opts(**self.cmbOpts).opts(**self.specOpts)

        #mobilogram
        rangeMz = streams.RangeY(y_range=(400, 1200), source=heatRt)
        if self.decimateMob:
            mobRt = decimate(hv.DynamicMap(self.plotMobilogramAcrossRT, streams=[rangeMz, framePointer]).opts(**self.cmbOpts).opts(**self.mobOpts))
        else:
            mobRt = hv.DynamicMap(self.plotChromatogram, streams=[rangeMz, framePointer]).opts(**self.cmbOpts).opts(**self.mobOpts)

        chrom = hv.DynamicMap(self.plotChromatogram, streams=[rangeMz, rangeMob]).opts(**self.cmbOpts).opts(**self.chromOpts)

        ## frame pointer graph on the chromatogram
        chromPointer = hv.DynamicMap(self.plotChromatogramPointer, streams=[framePointer]).opts(**self.cmbOpts).opts(**self.pointerChromOpts)


        chromTog = chrom * chromPointer


        # now that chromatogram is created, change the framePointer source
        framePointer.source = chrom

        

        return (heatRt + specRt + mobRt + chromTog ).cols(2)



# this class is for getting an MS1 plotter object
# gets all of the MS1 frames with a slider for retention time
# Usage:
# ms1Plot = MS1Plotter("< path to .d")
# ms1Plot.plot()
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
    def plot(self, **kwargs):

        # fetch MS1 frames list
        self.frames = self.fetchMS1Frames()

        return super().plot(**kwargs)


# This class is used to plot a swath window with a slider across retention time
# Main usage 
# swath = SwathPlotter("< path to .d folder")
# swath.plot(1,1) # for the first swath first multiplex 
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
    # frameRange = tuple of start and end (non inclusive) frames to plot
    def plotChrom(self, swathNum, multiplex, frameRange=None, **kwargs):
        # fetch frame list

        self.frames = self.fetchSwathFrames(swathNum)

        if frameRange is not None:
            self.frames = [ i for i in self.frames if i >= frameRange[0] and i < frameRange[1] ]

        meta = self.getSwathMetaInfo(swathNum, multiplex)


        print("Meta Info On Swath:\n{}".format(meta))


        # set the bounds appropriate for the swath
        self.mobOpts['ylim'] = (0, 20000)

        self.heatOpts['xlim'] = (meta['ImEnd'], meta['ImBegin'])
        self.specOpts['ylim'] = (0, 10000)


        # update the chromgram options based on supplied frames
        self.chromOpts['xlim'] = (min(self.frames)-1, max(self.frames) + 1)
        self.pointerChromOpts['xlim'] = self.chromOpts['xlim']

        return super().plotChrom(**kwargs)

    # swathNum is the group number (usually between 1-16)
    # multiplex should be value between 1 and 4)
    # frameRange = tuple of start and end (non inclusive) frames to plot
    def plot(self, swathNum, multiplex, frameRange=None, **kwargs):
        # fetch frame list

        self.frames = self.fetchSwathFrames(swathNum)

        if frameRange is not None:
            self.frames = [ i for i in self.frames if i >= frameRange[0] and i < frameRange[1] ]

        meta = self.getSwathMetaInfo(swathNum, multiplex)


        print("Meta Info On Swath:\n{}".format(meta))


        # set the bounds appropriate for the swath
        self.mobOpts['ylim'] = (0, 20000)

        self.heatOpts['xlim'] = (meta['ImEnd'], meta['ImBegin'])
        self.specOpts['ylim'] = (0, 10000)


        return super().plot(**kwargs)
