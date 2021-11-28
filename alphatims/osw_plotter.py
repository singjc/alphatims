# This library contains an additional object used for plotting OSW data

# The reasoning behind creating a different object is so the previous object is not dependent on the graphing depenendencies 

import alphatims.plotting
import alphatims.osw_parser
import param
import holoviews as hv
import bokeh.io

# plots transition data using OSW
class OSWPlotter:
    """
        Helper class to plot OSW data 
        Adopted from: https://stackoverflow.com/questions/58400769/panel-param-fileinput-widget-and-param-depends-interaction

        For more information on processing DIA data with OpenSwath please see http://openswath.org/en/latest/
    """
 
    def __init__(self, oswFile, timsTOF):
        self.oswFile = oswFile
        self.timsTOF = timsTOF


    def plotPrecursor(self, precursorId, rank=1, lineplot_kwargs=dict(x_axis_label='rt', remove_zeros=True), heatmap_kwargs=dict(), **kwargs):
        """
            Given a precursor Id, plot precursor XIC
            precursorId (int) --> id corresponding to osw file of precursor
            rank (int) --> which peak group rank to take, default is 1
        """
 
        self.oswFile = self.oswFile.subset_data_for_precursor(precursorId)
        precursorMetaData = self.oswFile.oswfile_data_current_precursor_subset[self.oswFile.oswfile_data_current_precursor_subset['peak_group_rank'] == 1.0].iloc[0]
        precursor=PeptideExtraction(mz = precursorMetaData['precursor_mz'], im = precursorMetaData['IM'], rt = precursorMetaData['RT'] * 60, peptideSequence=precursorMetaData['FullPeptideName'], charge=int(precursorMetaData['precursor_charge']))

        return precursor.inspect_peptide(self.timsTOF, lineplot_kwargs=lineplot_kwargs, heatmap_kwargs=heatmap_kwargs, **kwargs)



    def plotPrecursorAndTransition(self, precursorId, rank=1, lineplot_kwargs=dict(x_axis_label='rt', remove_zeros=True), heatmap_kwargs=dict(), **kwargs):
        """
            Given a precursor Id, plot the precursor XIC with its corresponding fragments
            precursorId (int) --> id corresponding to osw file of precursor
            rank (int) --> which peak group rank to take, default is 1
        """
        self.oswFile = self.oswFile.subset_data_for_precursor(precursorId)
        precursorMetaData = self.oswFile.oswfile_data_current_precursor_subset[self.oswFile.oswfile_data_current_precursor_subset['peak_group_rank'] == 1.0].iloc[0]

        # get transition info
        self.oswFile.fetchTransitionsFromPrecursor(precursorId)

 
        precursor=PeptideExtraction(mz = precursorMetaData['precursor_mz'], im = precursorMetaData['IM'], rt = precursorMetaData['RT'] * 60, peptideSequence=precursorMetaData['FullPeptideName'], charge=int(precursorMetaData['precursor_charge']), fragment_mzs=dict(self.oswFile.osw_data_fragment_ions.values))

        return precursor.inspect_peptide(self.timsTOF, lineplot_kwargs=lineplot_kwargs, heatmap_kwargs=heatmap_kwargs, **kwargs)


        self.oswFile.subset_data_for_precursor(precursorId)
        precursorMetaData = self.oswFile.subset_data_for_precursor[ oswFile.oswfile_data_current_precursor_subset['peak_group_rank'] == 1.0].iloc[0]

        # get transition info
        self.oswFile.fetchTransitionsFromPrecursor(precursorId)

        precursor=PeptideExtraction(mz = precursorMetaData['precursor_mz'], im = precursorMetaData['IM'], rt = precursorMetaData['RT'], peptide=['FullPeptideName'], fragment_ion=dict(self.oswFile.osw_data_fragment_ions.values))

        precursor.inspect_peptide()

    
    def plotPrecursorSelection(self):
        """
            Plot information in the precursor selection object
        """
        pass

    def plotPeptideSelection(self):
        """
            Plot information in the Peptide selection object
        """
        pass


class PeptideExtraction(param.Parameterized):
    """
        This class is a parameterized class with information on the peptide to plot
    """
    # mandatory peptide properties (must be set while initiating)
    peptideSequence = param.String(allow_None=False)
    mz = param.Number(doc='m/z of peptide', allow_None=False)
    im = param.Number(doc='im of peptide, (in 1/k0)', allow_None=False)
    rt = param.Number(doc='peptide retention time (seconds)', allow_None=False)


    #optional peptide properties
    charge = param.Integer(doc='Charge of the peptide')
    fragment_mzs = param.Dict(doc='dictionary of fragment mzs in the format {fragmentName:fragmentMz}', default={})


    # extraction preferences (optional, already have default values)
    rtExtraction = param.Number(doc='retention time extraction width (in seconds) on either side of the peptide', default=30, allow_None=False)
    imExtraction = param.Number(doc='IM extraction width (in 1/k0) on either side of the peptide', default=0.05, allow_None=False)
    ppm = param.Number(doc='ppm extraction around the m/z', default=50, allow_None=False)

    
    def __init__(self, **params):
        super().__init__(**params)


    # this functions is modified from the tutorial nodebook
    def inspect_peptide(self, dia_data, heatmap=False, save=False, lineplot_kwargs=dict(x_axis_label='rt', remove_zeros=True), heatmap_kwargs=dict()):
        rt_slice = slice(
            self.rt - self.rtExtraction,
            self.rt + self.rtExtraction 
        )
        im_slice = slice(
            self.im - self.imExtraction,
            self.im + self.imExtraction
        )
        precursor_mz_slice = slice(
            self.mz / (1 + self.ppm / 10**6),
            self.mz * (1 + self.ppm / 10**6)
        )
        precursor_indices = dia_data[
            rt_slice,
            im_slice,
            0, #index 0 means that the quadrupole is not used
            precursor_mz_slice,
            "raw"
        ]
        if heatmap:
            precursor_heatmap = alphatims.plotting.heatmap(
                dia_data.as_dataframe(precursor_indices),
                x_axis_label="rt",
                y_axis_label="mobility",
                title="precursor",
                width=250,
                height=250,
                **heatmap_kwargs
            )
            overlay = precursor_heatmap
        else:
            precursor_xic = alphatims.plotting.line_plot(
                dia_data,
                precursor_indices,
                width=900,
                label="precursor",
                **lineplot_kwargs
            )
            overlay = precursor_xic
        for fragment_name, m in self.fragment_mzs.items():
            fragment_mz_slice = slice(
                m / (1 + self.ppm / 10**6),
                m * (1 + self.ppm / 10**6)
            )
            fragment_indices = dia_data[
                rt_slice,
                im_slice,
                precursor_mz_slice,
                fragment_mz_slice,
                "raw"
            ]
            if len(fragment_indices) > 0:
                if heatmap:
                    fragment_heatmap = alphatims.plotting.heatmap(
                        dia_data.as_dataframe(fragment_indices),
                        x_axis_label="rt",
                        y_axis_label="mobility",
                        title=f"{fragment_name}: {m:.3f}",
                        width=250,
                        height=250,
                    )
                    overlay += fragment_heatmap
                else:
                    fragment_xic = alphatims.plotting.line_plot(
                        dia_data,
                        fragment_indices,
                        width=900,
                        label=fragment_name,
                        **lineplot_kwargs
                    )
                    overlay *= fragment_xic.opts(muted=True)
        if not heatmap:
            overlay.opts(hv.opts.Overlay(legend_position='bottom'))
            overlay.opts(hv.opts.Overlay(click_policy='mute'))
            overlay = overlay.opts(show_legend=True)
            if save:
                hv.save(overlay, f"{self.peptideSequence}_{self.charge}.html")
        return overlay.opts(
            title=f"{self.peptideSequence}_{self.charge}"
        )
