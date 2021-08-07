import sys
from pathlib import Path

from ase.io import read, write
from ase.io.formats import ioformats
from ase.calculators.abc import GetOutputsMixin
from ase.calculators.calculator import BaseCalculator


def read_stdout(args, createfile=None):
    """Run command in tempdir and return standard output.

    Helper function for getting version numbers of DFT codes.
    Most DFT codes don't implement a --version flag, so in order to
    determine the code version, we just run the code until it prints
    a version number."""
    import tempfile
    from subprocess import Popen, PIPE
    with tempfile.TemporaryDirectory() as directory:
        if createfile is not None:
            path = Path(directory) / createfile
            path.touch()
        proc = Popen(args,
                     stdout=PIPE,
                     stderr=PIPE,
                     stdin=PIPE,
                     cwd=directory,
                     encoding='ascii')
        stdout, _ = proc.communicate()
        # Exit code will be != 0 because there isn't an input file
    return stdout


class SingleFileReader:
    def __init__(self, fmt):
        self.fmt = fmt

    def __call__(self, path):
        output = read(path, format=self.fmt)
        cache = output.calc
        # XXX This is specially for things that
        # return Atoms + SinglePoint calculator
        # (Espresso)
        return dict(cache.properties())


class CalculatorTemplate:
    def __init__(self, name, implemented_properties,
                 input_file, output_file, input_format, reader):
        self.name = name
        self.implemented_properties = implemented_properties

        # Generalize: We need some kind of Writer and Reader
        # to handle multiple files at a time.
        self.input_file = input_file
        self.output_file = output_file
        self.input_format = input_format
        self.reader = reader

    def write_input(self, directory, atoms, parameters, properties):
        fmt = ioformats[self.input_format]
        if 'properties' in fmt.write.__code__.co_varnames:
            parameters = dict(parameters)
            parameters['properties'] = properties

        directory.mkdir(exist_ok=True, parents=True)

        path = directory / self.input_file
        write(path, atoms, format=fmt.name, **parameters)
        return path

    def execute(self, profile, directory):
        profile.run(directory,
                    self.input_file,
                    self.output_file)

    def read_results(self, directory):
        return self.reader(directory / self.output_file)

    def __repr__(self):
        return 'CalculatorTemplate({})'.format(vars(self))


def get_espresso_template():
    from ase.calculators.espresso import Espresso
    infile = 'espresso.pwi'
    outfile = 'espresso.pwo'
    return CalculatorTemplate(
        name='espresso',
        implemented_properties=Espresso.implemented_properties,
        input_file=infile,
        output_file=outfile,
        input_format='espresso-in',
        reader=SingleFileReader('espresso-out'))


class GenericFileIOCalculator(BaseCalculator, GetOutputsMixin):
    def __init__(self, template, profile, directory='.', parameters=None):
        self.template = template
        self.profile = profile

        # Maybe we should allow directory to be a factory, so
        # calculators e.g. produce new directories on demand.
        self.directory = Path(directory)

        if parameters is None:
            parameters = {}
        self.parameters = dict(parameters)

        self.atoms = None
        self.results = {}
        # XXX We are very naughty and do not call super constructor!

    def set(self, *args, **kwargs):
        raise RuntimeError('No setting parameters for now, please.  '
                           'Just create new calculators.')

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, self.template.name)

    def write_input(self, atoms, properties, system_changes):
        # XXX for socketio compatibility; remove later
        self.template.write_input(self.directory, atoms,
                                  self.parameters, properties)

    @property
    def implemented_properties(self):
        return self.template.implemented_properties

    @property
    def name(self):
        return self.template.name

    def calculate(self, atoms, properties, system_changes):
        self.atoms = atoms.copy()

        directory = self.directory

        self.template.write_input(directory, atoms, self.parameters,
                                  properties)
        self.template.execute(self.profile, directory)
        self.results = self.template.read_results(directory)
        # XXX Return something useful?

    def _outputmixin_get_results(self):
        return self.results
