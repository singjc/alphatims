import logging
import param
import panel as pn
import sqlite3
import pandas as pd

def check_sqlite_table(con, table):
        table_present = False
        c = con.cursor()
        c.execute('SELECT count(name) FROM sqlite_master WHERE type="table" AND name="%s"' % table)
        if c.fetchone()[0] == 1:
            table_present = True
        else:
            table_present = False
        c.fetchall()

        return(table_present)

class OSWFile(param.Parameterized):
    """
        Class to Process OSW File
        Adopted from: https://stackoverflow.com/questions/58400769/panel-param-fileinput-widget-and-param-depends-interaction

        For more information on processing DIA data with OpenSwath please see http://openswath.org/en/latest/
    """
    
    # File Selection paramater object
    oswfile = param.FileSelector(allow_None=False, doc="(Required) path to '.osw' file processed by OpenSwathWorkflow and PyProphet", precedence=1) 

    def __init__(self, **params):
        print("hello")
        super().__init__(**params)
        self.oswfile_data = None
        self.oswfile_data_current_peptide_subset = None
        self.oswfile_data_current_peptide_charge_subset = None
        self.oswfile_data_current_peptide_charge_peak_subset = None
        self.RT = None
        self.leftWidth = None
        self.rightWidth = None
        self.IM = None
        self.IM_leftWidth = None
        self.IM_rightWidth = None
        self.precursor_mz = None
        self.peak_group_rank = None
        self.Intensity = None
        self.ms2_m_score = None
        self.ipf_m_score = None

    @param.depends('oswfile', watch=True)  
    def process_file(self):
        if self.oswfile is not "":
            logging.info( f'INFO: Processing file - {self.oswfile}' )
            # Initiate connection to file
            con = sqlite3.connect(self.oswfile)
            # Check for necessary tables
            score_ms2_present = check_sqlite_table(con, "SCORE_MS2")
            score_ipf_present = check_sqlite_table(con, "SCORE_IPF")
            if score_ipf_present and not check_sqlite_table(con, "UNIMOD_CODENAME_MAPPING"):
                logging.error( f'INFO: Provided OSW file does not contain UNIMOD_CODENAME_MAPPING table for augmenting IPF scores. - {self.oswfile}' )
                score_ipf_present = False

            if score_ms2_present and not score_ipf_present:
                query = '''
                        SELECT 
                            PRECURSOR.ID AS precursor_id,
                            PEPTIDE.ID AS peptide_id,
                            PEPTIDE.MODIFIED_SEQUENCE as FullPeptideName,
                            PEPTIDE.UNMODIFIED_SEQUENCE as Sequence,
                            PRECURSOR.PRECURSOR_MZ as precursor_mz,
                            PRECURSOR.CHARGE as precursor_charge,
                            FEATURE.ID as feature_id,
                            FEATURE.EXP_RT/60 as RT,
                            FEATURE.LEFT_WIDTH/60 as leftWidth,
                            FEATURE.RIGHT_WIDTH/60 as rightWidth,
                            FEATURE.EXP_IM as IM,
                            FEATURE.EXP_IM+0.05 as IM_leftWidth,
                            FEATURE.EXP_IM-0.05 as IM_rightWidth,
                            FEATURE_MS2.AREA_INTENSITY as Intensity,
                            SCORE_MS2.RANK as peak_group_rank,
                            SCORE_MS2.QVALUE as ms2_m_score
                        FROM PRECURSOR
                        INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
                        INNER JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
                        INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
                        INNER JOIN PEPTIDE ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID
                        LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
                        WHERE PRECURSOR.DECOY = 0
                        AND SCORE_MS2.RANK IS NOT NULL -- For some reasons some features aren't scored?
                        '''
                logging.info( f'INFO: Reading peak group-level results.' )
                self.oswfile_data = pd.read_sql_query(query, con).sort_values(["Sequence", "FullPeptideName", "precursor_charge", "peak_group_rank"])
            elif score_ms2_present and score_ipf_present:
                query = '''
                        SELECT 
                            PRECURSOR.ID AS precursor_id,
                            PEPTIDE.ID AS peptide_id,
                            PEPTIDE.MODIFIED_SEQUENCE as FullPeptideName,
                            PEPTIDE.UNMODIFIED_SEQUENCE as Sequence,
                            PRECURSOR.PRECURSOR_MZ as precursor_mz,
                            PRECURSOR.CHARGE as precursor_charge,
                            FEATURE.ID as feature_id,
                            FEATURE.EXP_RT/60 as RT,
                            FEATURE.LEFT_WIDTH/60 as leftWidth,
                            FEATURE.RIGHT_WIDTH/60 as rightWidth,
                            FEATURE.EXP_IM as IM,
                            FEATURE.EXP_IM+0.05 as IM_leftWidth,
                            FEATURE.EXP_IM-0.05 as IM_rightWidth,
                            FEATURE_MS2.AREA_INTENSITY as Intensity,
                            SCORE_MS2.RANK as peak_group_rank,
                            SCORE_MS2.QVALUE as ms2_m_score,
                            SCORE_IPF.QVALUE as ipf_m_score
                        FROM PRECURSOR
                        INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID
                        INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
                        INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
                        INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
                        LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
                        LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
                        LEFT JOIN  (
                            SELECT 
                                SCORE_IPF.FEATURE_ID,
                                SCORE_IPF.PEPTIDE_ID,
                                UNIMOD_CODENAME_MAPPING.UNIMOD_ID,
                                SCORE_IPF.PRECURSOR_PEAKGROUP_PEP,
                                SCORE_IPF.QVALUE,
                                SCORE_IPF.PEP
                            FROM SCORE_IPF
                            INNER JOIN UNIMOD_CODENAME_MAPPING ON UNIMOD_CODENAME_MAPPING.CODENAME_ID = SCORE_IPF.PEPTIDE_ID
                        ) AS SCORE_IPF ON (SCORE_IPF.FEATURE_ID = FEATURE.ID AND SCORE_IPF.UNIMOD_ID = PEPTIDE.ID)
                        INNER JOIN PEPTIDE AS PEPTIDE_IPF ON SCORE_IPF.PEPTIDE_ID = PEPTIDE_IPF.ID
						INNER JOIN UNIMOD_CODENAME_MAPPING ON UNIMOD_CODENAME_MAPPING.UNIMOD_ID = PEPTIDE.ID
                        WHERE PRECURSOR.DECOY = 0
                        AND SCORE_MS2.RANK IS NOT NULL -- For some reasons some features aren't scored?
                        '''
                logging.info( f'INFO: Reading peak group-level IPF results.' )
                self.oswfile_data = pd.read_sql_query(query, con).sort_values(["Sequence", "FullPeptideName", "precursor_charge", "peak_group_rank"])
            else:
                logging.error( f'INFO: Provided OSW file does not seem to be scored. - {self.oswfile}' )
            # Close connection to file
            con.close()
        else:
            self.oswfile_data = None
            self.oswfile_data_current_peptide_subset = None
            self.oswfile_data_current_peptide_charge_subset = None
            self.oswfile_data_current_peptide_charge_peak_subset = None
            self.RT = None
            self.leftWidth = None
            self.rightWidth = None
            self.IM = None
            self.IM_leftWidth = None
            self.IM_rightWidth = None
            self.precursor_mz = None
            self.peak_group_rank = None
            self.Intensity = None
            self.ms2_m_score = None
            self.ipf_m_score = None
            logging.error( f'INFO: Cleared loaded OSW file...' )
        return self
    
    def subset_data_for_peptide(self, peptide):
        logging.info( f'INFO: Subsetting Identification Results for peptide {peptide}.' )
        self.oswfile_data_current_peptide_subset = self.oswfile_data[ self.oswfile_data.FullPeptideName == peptide ]
        return self
    
    def subset_data_for_charge(self, charge):
        logging.info( f'INFO: Subsetting Identification Results for precursor charge {charge}.' )
        self.oswfile_data_current_peptide_charge_subset = self.oswfile_data_current_peptide_subset[ self.oswfile_data_current_peptide_subset.precursor_charge == charge ]
        return self

    def subset_data_for_feature_id(self, feature_id):
        logging.info( f'INFO: Subsetting Identification Results for feature id {feature_id}.' )
        self.oswfile_data_current_peptide_charge_peak_subset = self.oswfile_data_current_peptide_charge_subset[ self.oswfile_data_current_peptide_charge_subset.feature_id == feature_id ]
        return self

    def get_rt_feature_data(self):
        if self.oswfile_data_current_peptide_charge_peak_subset is not None:
            self.precursor_mz = self.oswfile_data_current_peptide_charge_peak_subset.iloc[0]['precursor_mz']
            self.peak_group_rank = self.oswfile_data_current_peptide_charge_peak_subset.iloc[0]['peak_group_rank']
            self.RT=self.oswfile_data_current_peptide_charge_peak_subset.iloc[0]['RT']
            self.leftWidth=self.oswfile_data_current_peptide_charge_peak_subset.iloc[0]['leftWidth']
            self.rightWidth=self.oswfile_data_current_peptide_charge_peak_subset.iloc[0]['rightWidth']
            self.IM = self.oswfile_data_current_peptide_charge_peak_subset.iloc[0]['IM']
            self.IM_leftWidth = self.oswfile_data_current_peptide_charge_peak_subset.iloc[0]['IM_leftWidth']
            self.IM_rightWidth = self.oswfile_data_current_peptide_charge_peak_subset.iloc[0]['IM_rightWidth']
            self.Intensity = self.oswfile_data_current_peptide_charge_peak_subset.iloc[0]['Intensity']
            self.ms2_m_score = self.oswfile_data_current_peptide_charge_peak_subset.iloc[0]['ms2_m_score']
            if "ipf_m_score" in self.oswfile_data_current_peptide_charge_peak_subset.columns.tolist():
                self.ipf_m_score = self.oswfile_data_current_peptide_charge_peak_subset.iloc[0]['ipf_m_score']
        return self
