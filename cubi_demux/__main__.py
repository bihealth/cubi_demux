"""Main entry point of the ``cubi-demux`` application."""


class MissingConfiguration(Exception):
    """Raised on missing configuration."""


class DemuxApp:
    """Application class."""

    def __init__(self, config):
        #: Configuration from command line arguments.
        self.config = config
        # Check configuration.
        self.check_config()

    def check_config(self):
        """Check configuration."""
        if self.config.get('input_dir'):
            self.config['cubi_demux']['input_dir'] = self.config['input_dir']
            del self.config['input_dir']
        if not self.config['sample_sheet']:
            raise MissingConfiguration(
                'Configuration setting "sample_sheet" is missing from the '
                'YAML configuration file and/or --sample-sheet is missing '
                'from the command line.')
        if not self.config['cubi_demux']['input_dir']:
            raise MissingConfiguration(
                'Configuration setting "input_dir" is missing from the '
                'YAML configuration file and/or --input-dir is missing '
                'from the command line.')
