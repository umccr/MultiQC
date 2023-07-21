# Initialise the logger
import logging
import re
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


MIN_CNT_TO_SHOW_ON_PLOT = 5


class DragenFragmentLength(BaseMultiqcModule):
    def add_fragment_length_hist(self):
        data_by_sample = defaultdict(dict)

        for f in self.find_log_files("dragen/fragment_length_hist"):
            f_data_by_sample = parse_fragment_length_hist_file(f)
            for sample, data in f_data_by_sample.items():
                sample_clean = self.clean_s_name(sample, f)
                if sample_clean in data_by_sample:
                    log.debug(f'Duplicate sample name found! Overwriting: {sample_clean}')
                data_by_sample[sample_clean] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)

        # Check we have data before proceeding
        if not data_by_sample:
            return set()

        # Write data to file
        self.write_data_file(data_by_sample, "dragen_frag_len")

        smooth_points = 300
        self.add_section(
            name="Fragment length hist",
            anchor="dragen-fragment-length-histogram",
            description="""
            Distribution of estimated fragment lengths of mapped reads per sample.
            Only points supported by at least {} reads are shown to prevent long flat tail.
            The plot is also smoothed down to showing 300 points on the X axis to reduce noise.
            """.format(
                MIN_CNT_TO_SHOW_ON_PLOT
            ),
            plot=linegraph.plot(
                data_by_sample,
                {
                    "id": "dragen_fragment_length",
                    "title": "Dragen: Fragment length hist",
                    "ylab": "Number of reads",
                    "xlab": "Fragment length (bp)",
                    "ymin": 0,
                    "xmin": 0,
                    "tt_label": "<b>{point.x} bp</b>: {point.y} reads",
                    "smooth_points": smooth_points,
                },
            ),
        )
        return data_by_sample.keys()


def parse_fragment_length_hist_file(f):
    """
    T_SRR7890936_50pc.fragment_length_hist.csv

    #Sample: N_SRR7890889
    FragmentLength,Count
    36,1
    37,0
    38,0
    39,0
    40,0
    41,1
    ...
    39203,0
    39204,0
    39205,1
    #Sample: T_SRR7890936_50pc
    FragmentLength,Count
    53,2
    54,0
    ...
    39316,0
    39317,1
    """

    data_by_sample = defaultdict(dict)
    sample = None
    for line in f.get('f').splitlines():
        if line.startswith('#Sample'):
            sample = line.split(' ', maxsplit=1)[1]
            assert sample is not None
            assert len(sample) > 0
            assert sample not in data_by_sample
        elif line == 'FragmentLength,Count':
            continue
        else:
            fragment_length_str, count_str = line.split(',')
            fragment_length = int(fragment_length_str)
            count = int(count_str)
            # Skip low counts to in order to truncate uninformative long distribution tail
            if count < MIN_CNT_TO_SHOW_ON_PLOT:
                continue
            data_by_sample[sample][fragment_length] = count
    return data_by_sample
