# -*- coding: utf-8 -*-
"""The "actual" code, glued together by app and Snakefile.
"""

import collections
import csv
import textwrap

from ruamel import yaml

from .utils import listify


class Barcode:
    """Representation of a barcode used in sequencing."""

    def __init__(self, name, seq):
        """Initialize with values"""
        #: Barcode name to use
        self.name = name or None
        #: Barcode sequence to use
        self.seq = seq or None


class Library:
    """Representation of one library, used for the sample sheet generation
    """

    def __init__(self, name, reference, barcode, lanes, barcode2=None):
        """Initialize with the given parameters"""
        #: Name of the library, used in file name
        self.name = name
        #: UCSC ID for reference, used (and ignored) in v1 sample sheets
        self.reference = reference
        #: First barcode, usually used
        self.barcode = barcode
        #: List of integer lanes
        self.lanes = list(lanes)
        #: Second barcode, only used for dual indexing
        self.barcode2 = barcode2

    def barcode_seq(self, version):
        """Return barcode sequence for sample sheet"""
        assert version in [1, 2]
        if version == 1:
            if self.barcode2:
                return '{}-{}'.format(self.barcode.seq, self.barcode2.seq)
            else:
                return self.barcode.seq
        else:  # version == 2
            if self.barcode2:
                return (self.barcode.seq, self.barcode2.seq)
            else:
                return (self.barcode.seq,)

    def file_names(self, version, is_paired, lane=None, seq=None,
                   name=None):
        """Return list of file names that will be generated for this library

        ``version`` gives the RTA and demultiplexing program version (1 or 2)
        and ``is_paired`` gives the pairing flag
        """
        assert version in [1, 2]
        indices = [self.barcode.seq
                   if self.barcode and self.barcode.seq else 'NoIndex']
        reads = ('R1', 'R2') if is_paired else ('R1',)
        lanes = ['L{l:03d}'.format(l=l)
                 for l in self.lanes
                 if (lane is None) or (l == lane)]
        seq = seq
        if seq is None:
            seq = ''
        else:
            seq = seq + '_'
        if version == 1:
            tpl = '{sample_name}_{index}_{lane}_{read}_001.fastq.gz'
        else:
            tpl = '{sample_name}_{seq}{lane}_{read}_001.fastq.gz'
        return list(sorted([
            tpl.format(sample_name=name or self.name, index=index, lane=lane,
                       read=read, seq=seq)
            for index in indices
            for read in reads
            for lane in lanes]))


NameTokens = collections.namedtuple(
    'NameTokens',
    ['date', 'instrument', 'run', 'slot', 'flow_cell', 'label'])


class FlowCell:
    """Representation of one flow cell, used for the sample sheet generation
    """

    def __init__(self, name, num_lanes, operator, rta_version, is_paired,
                 read_length, libraries=[]):
        #: Name of the flow cell
        self.name = name
        #: Number of lanes
        self.num_lanes = num_lanes
        #: Tokens from the flow cell
        self.name_tokens = self._build_name_tokens(name)
        #: Name of the operator
        self.operator = operator
        #: RTA version, currently only 1 and 2 are accepted
        assert rta_version in {1, 2}
        self.rta_version = rta_version
        #: Flag indicating whether the run was done in paired-end moe
        self.is_paired = is_paired
        #: Read length
        self.read_length = read_length
        #: ``Library`` objects
        self.libraries = list(libraries)

    @property
    def undetermined_libraries(self):
        """Return Library objects for undetermined libraries of lanes
        that we have data for"""
        lanes = set()
        for lib in self.libraries:
            lanes |= set(lib.lanes)
        result = []
        for lane in lanes:
            result.append(Library(
                'lane{}'.format(lane), 'none',
                Barcode('Undetermined', 'Undetermined'), [lane]))
        return result

    def _build_name_tokens(self, name):
        """Build Flow Cell name tokens"""
        tokens = name.split('_')
        date = tokens[0]
        instrument = tokens[1]
        run = tokens[2]
        if len(tokens[3]) == 1:
            slot = tokens[3]
            flow_cell = tokens[4]
            label = '_'.join(tokens[5:])
        else:
            slot = None
            flow_cell = tokens[3]
            label = '_'.join(tokens[4:])
        return NameTokens(date, instrument, run, slot, flow_cell, label)


class SampleSheet:
    """Class for loading ``FlowCell`` from YAML and writing out as sample
    sheet
    """

    @classmethod
    def from_yaml_data(klass, data):
        """Load from configuration 'sample_sheet' entry and return
        ``SampleSheet``.
        """
        flow_cell = FlowCell(
            data['name'], data['num_lanes'], data['operator'],
            data['rta_version'], data['is_paired'], data['read_length'])
        for d in data['libraries']:
            barcode = Barcode(d['barcode']['name'], d['barcode']['seq'])
            if 'barcode2' in d:
                barcode2 = Barcode(d['barcode2']['name'], d['barcode2']['seq'])
            else:
                barcode2 = None
            flow_cell.libraries.append(Library(
                d['name'], d['reference'], barcode, d['lanes'], barcode2))
        return SampleSheet(flow_cell)

    def __init__(self, flow_cell):
        """Initialize the sample sheet with the given ``FlowCell`` object"""
        self.flow_cell = flow_cell

    def print_v1(self, f):
        """Write sample sheet to file-like object f in bcl2fastq v1 format"""
        writer = csv.writer(f)
        header = ['FCID', 'Lane', 'SampleID', 'SampleRef', 'Index',
                  'Description', 'Control', 'Recipe', 'Operator',
                  'SampleProject']
        writer.writerow(header)
        recipe = 'PE_indexing' if self.flow_cell.is_paired else 'SE_indexing'
        for lib in self.flow_cell.libraries:
            for lane in sorted(lib.lanes):
                data = [self.flow_cell.name_tokens.flow_cell, str(lane),
                        lib.name, lib.reference, lib.barcode_seq(1), '', 'N',
                        recipe, self.flow_cell.operator, 'Project']
                writer.writerow(data)

    def print_v2(self, f):
        """Write sample sheet to file-like object f in bcl2fastq v2 format"""
        writer = csv.writer(f)
        # Write [Data] Section
        writer.writerow(['[Data]'])
        barcode2 = any([library.barcode2
                        for library in self.flow_cell.libraries])
        if barcode2:
            writer.writerow(['lane', 'sample_id', 'index', 'index2',
                             'sample_project'])
        else:
            writer.writerow(['lane', 'sample_id', 'index', 'sample_project'])
        rows = []
        for lib in self.flow_cell.libraries:
            for lane in sorted(lib.lanes):
                rows.append([lane, lib.name] + list(lib.barcode_seq(2)) +
                            ['Project'])
        for row in sorted(rows):
            writer.writerow(list(map(str, row)))


class Bcl2FastqBaseWrapper:
    """Base class for Bcl2FastqVXWrapper classes"""

    def __init__(self, config):
        #: Configuration.
        self.config = config
        #: Construct SampleSheet object from data previously loaded from YAML.
        self.sample_sheet = SampleSheet.from_yaml_data(
            config['sample_sheet'][0])
        #: FlowCell instance to use.
        self.flow_cell = self.sample_sheet.flow_cell

    def barcode_mismatches(self):
        """Return barcode mismatches (configuration or version-specific
        defaults).
        """
        raise NotImplementedError()

    def get_tiles_arg(self):
        """Return --tiles arg if any"""
        demux_config = self.config['cubi_demux']
        # Check whether all lanes are specified in the sheet
        lanes = set()
        for lib in self.flow_cell.libraries:
            lanes |= set(lib.lanes)
        all_lanes_in_sheet = (len(lanes) == self.flow_cell.num_lanes)
        # Shortcut to whether any args are given for tiles
        args_given = bool(demux_config['tiles'] or demux_config['lanes'])
        # Handle case of no args given and fake --lane parameter if not all
        # lanes mentioned in sheet.
        if not args_given and all_lanes_in_sheet:
            return ''  # no args
        elif not args_given and not all_lanes_in_sheet:
            demux_config['lanes'] = list(lanes)
        # Handle case of selected lanes or selected tiles
        if demux_config['tiles']:
            return '--tiles {}'.format(','.join(demux_config['tiles']))
        else:
            regexes = ['s_{}'.format(lane)
                       for lane in sorted(demux_config['lanes'])]
            return '--tiles {}'.format(','.join(regexes))

    def lane_enabled(self, lane):
        """Return whether lane enabled via args"""
        if not self.config['cubi_demux']['lanes']:
            return True
        else:
            return lane in self.config['cubi_demux']['lanes']


class Bcl2FastqV1Wrapper(Bcl2FastqBaseWrapper):
    """Implementation of demultiplexing for bcl2fastq 1.x"""

    def barcode_mismatches(self):
        """Return barcode mismatches (configuration or V1 defaults)."""
        if self.config['cubi_demux']['barcode_mismatches'] is None:
            return 1
        else:
            return self.config['cubi_demux']['barcode_mismatches']

    @listify
    def get_output_files(self):
        """Return list of paths to output files"""
        flow_cell = self.parent.sample_sheet.flow_cell
        for lib in flow_cell.libraries + flow_cell.undetermined_libraries:
            for lane in sorted(lib.lanes):
                if not self.lane_enabled(lane):
                    continue  # skip
                sample_name = lib.name
                if lib.barcode_seq(1) == 'Undetermined':
                    sample_name = 'Undetermined'
                out_dir = ('{args.output_dir}/{sample_name}/'
                           '{flowcell}/L{lane:03d}').format(
                               flowcell=flow_cell.name_tokens.flow_cell,
                               args=self.parent.args,
                               sample_name=sample_name,
                               lane=lane)
                for fname in lib.file_names(flow_cell.rta_version,
                                            flow_cell.is_paired,
                                            lane):
                    yield '{out_dir}/{fname}'.format(out_dir=out_dir,
                                                     fname=fname)

    def print_sheet(self, f):
        """Print sheet to file ``f``."""
        self.sample_sheet.print_v1(f)


class Bcl2FastqV2Wrapper(Bcl2FastqBaseWrapper):
    """Implementation of demultiplexing for bcl2fastq 2.x"""

    def __init__(self, config):
        super().__init__(config)
        #: Sample mapping from library name to SXX id.
        self.sample_map = self._build_sample_map()

    def _build_sample_map(self):
        """Build mapping from library name to SXX id"""
        result = {}
        rows = [
            (lane, lib.name)
            for lib in self.flow_cell.libraries
            for lane in lib.lanes]
        i = 1
        for _, name in sorted(set(rows)):
            if name not in result:
                result[name] = 'S{}'.format(i)
                i += 1
        return result

    def barcode_mismatches(self):
        """Return barcode mismatches (configuration or V2 defaults)."""
        if self.config['cubi_demux']['barcode_mismatches'] is None:
            return 0
        else:
            return self.config['cubi_demux']['barcode_mismatches']

    @listify
    def get_output_files(self):
        """Return list of paths to output files"""
        flow_cell = self.sample_sheet.flow_cell
        for lib in flow_cell.libraries + flow_cell.undetermined_libraries:
            for lane in sorted(lib.lanes):
                if not self.lane_enabled(lane):
                    continue  # skip
                sample_name = lib.name
                if lib.barcode_seq(1) == 'Undetermined':
                    sample_name = 'Undetermined'
                output_dir = self.config['cubi_demux']['output_dir']
                out_dir = ('{output_dir}/{sample_name}/'
                           '{flowcell}/L{lane:03d}').format(
                               flowcell=flow_cell.name_tokens.flow_cell,
                               output_dir=output_dir,
                               sample_name=sample_name,
                               lane=lane)
                seq = self.sample_map.get(sample_name, 'S0')
                name = ('Undetermined' if lib.barcode_seq(1) == 'Undetermined'
                        else lib.name)
                for fname in lib.file_names(flow_cell.rta_version,
                                            flow_cell.is_paired,
                                            lane, seq, name):
                    yield '{out_dir}/{fname}'.format(out_dir=out_dir,
                                                     fname=fname)

    def print_sheet(self, f):
        """Print sheet to file ``f``."""
        self.sample_sheet.print_v2(f)
