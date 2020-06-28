

from LFPyStim import TemplateCell
import aberraAxon

import neuron


class MyelinatedTemplateCell(TemplateCell):

    def __init__(self, templatefile='LFPyCellTemplate.hoc',
                 templatename='LFPyCellTemplate',
                 templateargs=None,
                 verbose=True,
                 **kwargs):

        TemplateCell.__init__(self,
                 templatefile=templatefile,
                 templatename=templatename,
                 templateargs=templateargs  ,
                 verbose=verbose,
                 **kwargs)

    def add_axons(self):
        myelinater = aberraAxon.MyelinatedCell(hocObj = neuron.h)
        

        myelinater.loadcell(self.templatename, myelinate_ax=True,  synapses = False)
        
        self.template = myelinater.cell

        for sec in neuron.h.allsec():
            if neuron.h.ismembrane("xtra", sec=sec):
                for x in sec:
                    x.es_xtra = 0


        neuron.h._ref_stim_xtra[0] = 0




    def _load_geometry(self):
        """Load the morphology-file in NEURON"""
        """ This is a slightly modified version from parent class, that modifies axon morphology"""

        try:
            neuron.h.sec_counted = 0
        except LookupError:
            neuron.h('sec_counted = 0')

        # the python cell object we are loading the morphology into:
       
        self.add_axons()

        # perform a test if the morphology is already loaded:
        seccount = 0
        for sec in self.template.all:
            seccount += 1
        if seccount == 0:
            # import the morphology, try and determine format
            fileEnding = self.morphology.split('.')[-1]

            if not fileEnding == 'hoc' or fileEnding == 'HOC':
                # create objects for importing morphologies of different formats
                if fileEnding == 'asc' or fileEnding == 'ASC':
                    Import = neuron.h.Import3d_Neurolucida3()
                    if not self.verbose:
                        Import.quiet = 1
                elif fileEnding == 'swc' or fileEnding == 'SWC':
                    Import = neuron.h.Import3d_SWC_read()
                elif fileEnding == 'xml' or fileEnding == 'XML':
                    Import = neuron.h.Import3d_MorphML()
                else:
                    raise ValueError('%s is not a recognised morphology file format! ').with_traceback('Should be either .hoc, .asc, .swc, .xml!'
                                                                                                       % self.morphology)

                # assuming now that morphology file is the correct format
                try:
                    Import.input(self.morphology)
                except:
                    if not hasattr(neuron, 'neuroml'):
                        raise Exception('Can not import, try and copy the ' +
                                        'nrn/share/lib/python/neuron/neuroml ' +
                                        'folder into %s' % neuron.__path__[0])
                    else:
                        raise Exception(
                            'something wrong with file, see output')
                try:
                    imprt = neuron.h.Import3d_GUI(Import, 0)
                except:
                    raise Exception('See output, try to correct the file')

                # instantiate the cell object
                if fileEnding == 'xml' or fileEnding == 'XML':
                    # can not currently assign xml to cell template
                    try:
                        imprt.instantiate(self.template)
                    except:
                        raise Exception("this xml file is not supported")
                else:
                    imprt.instantiate(self.template)

            else:
                neuron.h.execute("xopen(\"%s\")" %
                                 self.morphology, self.template)

        # set shapes and create sectionlists
    
        print("added_axons")
        neuron.h.define_shape()
        print("shape defined")
        self._create_sectionlists()
