"""
This is a re-implementation of the axon models preseneted in:

Aberra AS, Peterchev AV, Grill WM (2018) Biophysically realistic neuron models 
for simulation of cortical stimulation. J Neural Eng 15:066023 

This a near exact python port of the original hoc code found here:
https://senselab.med.yale.edu/modeldb/ShowModel?model=241165

Reasons to make this that I am not a hoc expert and I wanted access to these models in a way that 
allowed me to be able to modify them flexibly in order to integrate with other libraries.

As an exercise porting this algorithm was also valuable in learning how to use the python-hoc NEURON interface.


Huge thanks to the authors of the original code and methods as all ideas here are theirs.
This is mostly just a change from hoc syntax to python syntax for many functions


Author of original source code: Aman Aberra

Author of python port: Liam Long
For questions specific to port Contact: liam.long@mail.utoronto.ca

For questions on the neuroscience or algorithm please consult paper cited above
"""


import os
import pdb
import sys
import shutil

import neuron


import numpy as np

import matplotlib
matplotlib.use('TkAgg')

class HocUndefinedError(Exception):
    pass

class strdef():

    debug = True

    def __init__(self, var, hocObj=None):


        assert(type(var) == str)
        self.name = var

        if hocObj is None:
            raise HocUndefinedError
            self.hocObj = neuron.h
        else:
            self.hocObj = hocObj

        # initialize objref in hoc
        try:
            self.create(var)

        except Exception as e:
            print(e)
            if self.debug:
                pdb.set_trace()

    def set(self, toWhat):
        

        self.hocObj("{} = \"{}\"".format(self.name, toWhat))
        try:
            assert(getattr(self.hocObj, self.name) == toWhat)
        except AssertionError as e:
            print(e)
            pdb.set_trace()

    def get(self):
        return getattr(self.hocObj, self.name)

    def create(self, var):
        self.hocObj("strdef %s" % var)


def create(var, size=None, hocObj=None):
    assert(type(var) == str)
    assert(type(size) == str)


    if hocObj is None:
        raise HocUndefinedError
        hocObj = neuron.h

    if size is None:
        hocObj("create {}".format(var))
    else:
        hocObj("create {}[{}]".format(var, size))
    return getattr(hocObj, var)


def refpointer(var, val, hocObj=None):
    assert(type(var) == str)
    if hocObj is None:
        raise HocUndefinedError
        hocObj = neuron.h
    hocObj("{} = {}" .format(var, val))
    return getattr(hocObj, "_ref_{}".format(var))


def objref(var, hocObj=None):
    if hocObj is None:
        raise HocUndefinedError
        hocObj = neuron.h

    assert(type(var) == str)
    hocObj("objref {}".format(var))
    return getattr(hocObj, var)


class interpCoordinates():

    debug = False

    def __init__(self, caller):
        self.caller = caller

        self.xx = objref("xx", hocObj = self.caller.hocObj)
        self.yy = objref("yy", hocObj = self.caller.hocObj)
        self.zz = objref("zz", hocObj = self.caller.hocObj)
        self.length = objref("length", hocObj = self.caller.hocObj)

        self.xint = objref("xint", hocObj = self.caller.hocObj)
        self.yint = objref("yint", hocObj = self.caller.hocObj)
        self.zint = objref("zint", hocObj = self.caller.hocObj)
        self.range = objref("range", hocObj = self.caller.hocObj)

        self.current_dir = strdef("current_dir", hocObj = self.caller.hocObj)
        self.current_dir.set(self.caller.hocObj.getcwd())

        self.xList = self.caller.hocObj.List()
        self.yList = self.caller.hocObj.List()
        self.zList = self.caller.hocObj.List()

        self.numSect = refpointer("numSect", "1", hocObj = self.caller.hocObj)
        self.numComp = refpointer("numSect", "1", hocObj = self.caller.hocObj)

        self.secrefs = objref("secref", hocObj = self.caller.hocObj)

    def getcoords(self):
        print("get coords.....")
        self.xList.remove_all()
        self.yList.remove_all()
        self.zList.remove_all()
        
        for sec in self.caller.cell.all:
         

            if self.caller.hocObj.ismembrane("xtra"):

                nn = int(self.caller.hocObj.n3d())
                
                self.xx = self.caller.hocObj.Vector(nn)
                self.yy = self.caller.hocObj.Vector(nn)
                self.zz = self.caller.hocObj.Vector(nn)
                self.length = self.caller.hocObj.Vector(nn)

                for ii in range(nn):
                    self.xx.x[ii] = self.caller.hocObj.x3d(ii)
                    self.yy.x[ii] = self.caller.hocObj.y3d(ii)
                    self.zz.x[ii] = self.caller.hocObj.z3d(ii)
                    self.length.x[ii] = self.caller.hocObj.arc3d(ii)
             
                self.length.div(self.length.x[nn - 1])

                self.range = self.caller.hocObj.Vector(self.caller.hocObj.nseg + 2)
                self.range.indgen(1 / self.caller.hocObj.nseg)
                self.range.sub(1 / (2 * self.caller.hocObj.nseg))
                self.range.x[0] = 0
                self.range.x[self.caller.hocObj.nseg + 1] = 1

                self.xint = self.caller.hocObj.Vector(self.caller.hocObj.nseg + 2)
                self.yint = self.caller.hocObj.Vector(self.caller.hocObj.nseg + 2)
                self.zint = self.caller.hocObj.Vector(self.caller.hocObj.nseg + 2)

                self.xint.interpolate(self.range, self.length, self.xx)
                self.yint.interpolate(self.range, self.length, self.yy)
                self.zint.interpolate(self.range, self.length, self.zz)

                # original code x_xtra(xr) = .....
                # x_xtra is xtra mechanism of cas
                # the below should give pointer to that mechanism at same location
                for ii in range(1, self.caller.hocObj.nseg + 1):
                    xr = self.range.x[ii]
                    sec(xr).x_xtra = self.xint.x[ii]
                    sec(xr).y_xtra = self.yint.x[ii]
                    sec(xr).z_xtra = self.zint.x[ii]

                self.xint.remove(self.caller.hocObj.nseg + 1)
                self.xint.remove(0)
                self.yint.remove(self.caller.hocObj.nseg + 1)
                self.yint.remove(0)
                self.zint.remove(self.caller.hocObj.nseg + 1)
                self.zint.remove(0)

                self.xList.append(self.xint)
                self.yList.append(self.yint)
                self.zList.append(self.zint)

        self.numSect = self.xList.count()
        print("numSect was %i\n" % self.numSect)

    def getSecRefs(self):
        self.getcoords()
        self.secrefs = self.caller.hocObj.List()
        secnum = 0
        self.numComp = 0

        for sec in self.caller.cell.all:
            if self.caller.hocObj.ismembrane("xtra"):
                self.secrefs.append(self.caller.hocObj.SectionRef())  # SectionRef?????????
                secnum += 1
                self.numComp = self.numComp + self.caller.hocObj.nseg
        print("Created List of SectionRefs for %g sections\n" % secnum)

        self.assign_section_types()
        self.assign_order2()


    def traverse_tree(self, order, current_secref):

        current_secref.sec.order_xtra = order

        if current_secref.nchild() == 0:
            if self.debug:
                print("Reached terminal at {}, order is {}".format(self.caller.hocObj.secname(sec = current_secref.sec), 
                    current_secref.sec.order_xtra))

            #return to branch point
            not_branch = 1
            while not_branch:
               
                current_secref = self.caller.hocObj.SectionRef(sec = current_secref.parent)
                if current_secref.nchild() > 1:
                    not_branch = 0
                   
                else:
                    order = order - 1
               
        elif current_secref.nchild() == 1:
            if self.debug:
                print("1 Children in {}, order is {}".format(current_secref.nchild(), self.caller.hocObj.secname(sec = current_secref.sec), 
                    current_secref.sec.order_xtra))
            
          
            current_secref = self.caller.hocObj.SectionRef(sec = current_secref.child[0] )
       

            self.traverse_tree(order, current_secref)

        else:
            if self.debug:
                rint("{} Children in {}, order is {}".format(current_secref.nchild(),
                self.caller.hocObj.secname(sec=current_secref.sec),current_secref.sec.order_xtra) )


            order = order + 1
            for j in range(current_secref.nchild()): 
           
            

                self.traverse_tree(order, self.caller.hocObj.SectionRef(sec = current_secref.child[j]))

 

    def assign_order2(self):
        #  self.caller.hocObj.object_push(self.secrefs[0].sec)
        i = 1
        for n in range(self.secrefs[0].nchild()):

            print("children of soma:    %i"%n)
            oseci_secref = self.caller.hocObj.SectionRef(sec = self.secrefs[0].child[n])
            self.traverse_tree(i, oseci_secref)
     

    def assign_section_types(self):
        print("numsect: ",self.numSect)
        for i in range(self.numSect): 
            if self.secrefs[i].has_trueparent() == 0:
              
                self.secrefs[i].sec.type_xtra = 1
            else:
               
              

                parent_sec = self.caller.hocObj.SectionRef(sec = self.secrefs[i].parent().sec)
                parent_nchildren = parent_sec.nchild()

                Li = self.caller.hocObj.L

                if self.caller.hocObj.ismembrane("pas", sec=self.secrefs[i].sec):

                    Lambdai = self.Lambda(self.caller.hocObj.cas(), 1)*1000
                else:
                    Lambdai = self.caller.hocObj.L

                nn = self.caller.hocObj.nseg


                if parent_nchildren == 1:  # intermediate
                    if self.secrefs[i].nchild() == 0: # termination
                        for ix in range(1, nn+1):
                            self.secrefs[i].sec((2*ix-1)/(2*nn)).type_xtra = 3
                        self.secrefs[i].sec(1).type_xtra = 2

                    elif self.secrefs[i].nchild() == 1: # intermediate
                        for ix in range(1, nn+1):
                            self.secrefs[i].sec((2*ix-1)/(2*nn)).type_xtra = 3
                        
                    else:
                        for ix in range(1, nn+1): # bifurcate
                            self.secrefs[i].sec((2*ix-1)/(2*nn)).type_xtra = 3
                        self.secrefs[i].sec(1).type_xtra = 4

                else:  # bifurcation
                    if self.secrefs[i].nchild() == 0:
                        if Li <= Lambdai:
                            for ix in range(1, nn+1):
                                self.secrefs[i].sec((2*ix-1)/(2*nn)).type_xtra = 3
                            self.secrefs[i].sec(0).type_xtra = 6
                            self.secrefs[i].sec(1).type_xtra = 5
                        else:
                            for ix in range(1, nn+1):
                                self.secrefs[i].sec((2*ix-1)/(2*nn)).type_xtra = 3
                            self.secrefs[i].sec(0).type_xtra = 6
                            self.secrefs[i].sec(1).type_xpatra = 2

                    elif self.secrefs[i].nchild() == 1:
                        for ix in range(1, nn+1):
                            self.secrefs[i].sec((2*ix-1)/(2*nn)).type_xtra = 3
                        self.secrefs[i].sec(0).type_xtra = 6
                    else:
                        if Li <=Lambdai:
                            for ix in range(1,nn+1):
                                self.secrefs[i].sec((2*ix-1)/(2*nn)).type_xtra = 3
                            self.secrefs[i].sec(0).type_xtra = 6
                            self.secrefs[i].sec(1).type_xtra = 7
                        else:
                            for ix in range(1, nn+1):
                                self.secrefs[i].sec((2 * ix - 1)/(2 * nn)).type_xtra = 3
                            self.secrefs[i].sec(0).type_xtra = 6
                            self.secrefs[i].sec(1).type_xtra = 4
        print("Assigned section types to each section in type_xtra")






        del parent_sec

    def Lambda(self, sec, pos): 

        l = 10 * np.sqrt(((1 / sec(pos).g_pas) * sec.diam * 1e-4) / (4 * sec.Ra))

        return l


class CellLoader():

    debug = True

    def __init__(self, caller):
        self.caller = caller
        self.cell_id = refpointer("cell_id", "0", hocObj = self.caller.hocObj)

      
        self.cells = objref("cell", hocObj = self.caller.hocObj)
        self.nil = objref("nil", hocObj = self.caller.hocObj)

        self.main_ax_list = objref("main_ax_list", hocObj = self.caller.hocObj)
        self.strobj = objref("strobj", hocObj = self.caller.hocObj)

        self.secnames = objref("secnames", hocObj = self.caller.hocObj)

        self.current_dir = strdef("current_dir" , hocObj = self.caller.hocObj)
        self.current_dir.set("getcwd()")

        self.tstr2 = strdef("tstr2", hocObj = self.caller.hocObj)
        self.cell_dir = strdef("cell_dir", hocObj = self.caller.hocObj)
        self.createsim_file = strdef("createsim_file", hocObj = self.caller.hocObj)


    def cell_chooser(self, cell_name, myelinate_ax = True, load_synapses=True, loadedtemplate = None):
        self.cell_name = 'neurons/' + cell_name
        
      
        if loadedtemplate is None:
            #self.caller.hocObj("forall delete_section()")
            self.get_cell(load_synapses)
        else:
            self.cell = loadedtemplate
        print("+++++++++++++++++++++++++")
        print(self.cell)
        if myelinate_ax:
            self.add_axons(cell_name)
        else:
            self.setupExtracellStim() # make sure we can apply extrcellular stim
                                    
                

        print("cell loaded")


    def add_axons(self, cell_name, cell_id = 16):
        #TODO: fix this cell_id parameter
        print("+++++++++++++++++++++++++")
        print(self.cell)
        for sec in self.caller.cell.all:  
             # print("Inserting xtra into {}".format(self.caller.hocObj.secname(sec=sec)))

            sec.insert('xtra')
            sec.insert('extracellular')


        # Scaling
        # scale_axon_diameter = 1.322 # TODO figure how to set this paramter for now default in Aberra code
        # self.caller.Morphology.scale_diam2(scale_axon_diameter, self.caller.cellChooser.cell.axonal)


        self.caller.interpCoordinates.getSecRefs()

        meth1cells = self.caller.hocObj.Vector()
        meth1cells.append(6, 8, 20, 21, 22, 25)
        
        # TODO find a way to streamline this without cellids? 
        # for now defaults to ax2 as used for most cells in AberraetAl2018 original code
        if not meth1cells.contains(cell_id):
            self.main_ax_list = self.get_main_ax2()
            print("got main_ax 2")
        else:
            self.main_ax_list = self.get_main_ax()



        self.caller.Morphology.myelinate_axon(self.caller.CellLoader.cell.axonal)
      

        self.setupExtracellStim()
    
        if not meth1cells.contains(cell_id):
            self.main_ax_list = self.get_main_ax2()
            print("got main_ax 2")
        else:
            self.main_ax_list = self.get_main_ax()
        print("defining shape")
        self.caller.hocObj.define_shape()
        print("returning from axons")

    def setupExtracellStim(self):
        numComp = 0 # numComp not needed?
        for sec in self.caller.cell.all:
           
            sec.insert("xtra")
            sec.insert("extracellular")
            numComp = sec.nseg + numComp
        

        print("Inserted xtra and extracellular in all {} compartments\n".format(numComp))
        self.setpointers()


    def setpointers(self):
        self.caller.interpCoordinates.getSecRefs() 
        print("setting pointers?")
        for sec in self.caller.cell.all:


            if self.caller.hocObj.ismembrane("xtra", sec=sec) and self.caller.hocObj.ismembrane("extracellular", sec=sec):
                # should be a way to do this without resortin to this ugly str coding but fine for now
                self.caller.hocObj("for (x,0){ setpointer ex_xtra(x), e_extracellular(x) }")

    
    

    def getTemplate(self, directory=''):
        try: 
            f = open(directory + 'template.hoc')
        except:
            f = open('template.hoc')
        for line in f:
            if 'begintemplate' in line:
                templatename = line.split(' ')[1][:-1]
                # print(templatename)
                return templatename

    def replace_line(self, original, newname, origline, newline):
        print(original, newname, origline, newline)
        shutil.copy(original, newname)
        if not isinstance(origline, list):
            origline = [origline]
        if not isinstance(newline, list):
            newline = [newline]
        with open(newname, mode='w') as new_f:
            with open(original, mode='r') as old_f:
                for line in old_f.readlines():
                    for orig, new in zip(origline,newline):
                        line = line.replace(orig, new)
                    new_f.write(line)
      

    def gen_template_withaxon(self, template, synapses = False):
        axontemplate = 'axon_' + template
        
        axon_removal_string = '    replace_axon()'
        self.replace_line(template, axontemplate, axon_removal_string, '//' + axon_removal_string)

        removal_string = 'forall delete_section()'
        self.replace_line(template, axontemplate, removal_string, '//' + removal_string)
   
        return axontemplate

    def gen_synapses(self, templatefile, directory = ''):

        """
        if using synapses we need modify how they attach to axons

        For now let's say that they attach to axon segment 0 (which does not get deleted)
        """
        newtemplate = "syn_" + templatefile
        self.replace_line(templatefile, newtemplate, r'load_file("synapses/synapses.hoc")', r'load_file("synapses/axon_synapses.hoc")')
        
        template = 'synapses/synapses.hoc'
        new = 'synapses/axon_synapses.hoc'
        synload_string = r'synapse_file = new File("synapses/synapses.tsv")'
        new_string =  r'synapse_file = new File("synapses/axon_synapses.tsv")'
        print(synload_string, new_string)
     
        
        synload_string2 = r'printf'
        new_string2 =  r'//printf'
        synload_string3 = r' synapse_type_name, synapse_id, pre_cell_id, id_mtype_map.o(pre_mtype).s, sectionlist_name'
        new_string3 =  r'// synapse_type_name, synapse_id, pre_cell_id, id_mtype_map.o(pre_mtype).s, sectionlist_name'
        synload_string4 = r'sectionlist_index, seg_x, dep'
        new_string4 =  r'//sectionlist_index, seg_x, dep'

        old_string =  [synload_string]
        new_string =  [new_string]

        print(synload_string2, new_string2)
        self.replace_line(template, new, old_string, new_string)


        template = 'synapses/synapses.tsv'
        new = 'synapses/axon_synapses.tsv'
        print(synload_string, new_string)
        with open(directory + template, mode = 'r') as old_f:
            with open(directory + new, mode = 'w') as new_f:
                for line in old_f.readlines():
                    params = line.split()
                    if len(params) != 13:
                        new_f.write(line)
                        continue
                    isAxon = params[3] == '3'
                    if isAxon:
                       
                        params[4] = '0' # set section number to axon[0]
                        newline= ''
                        for param in params:
                           
                            newline += param+'\t'
                      
                        line = newline[:-1] + '\n' # replace last space with newline char

                    new_f.write(line)

        return newtemplate

    def get_cell(self, synapses):
        """Create the cell model"""
        self.cell_dir.set(self.cell_name)
        try:
            os.chdir(self.cell_dir.get())
        except:
            pass
        templatename = self.getTemplate(directory='')
        templatefile = self.gen_template_withaxon('template.hoc')

        if synapses:
            templatefile = self.gen_synapses(templatefile)

    
        # Load main cell template
        print(templatefile)
        self.caller.hocObj.load_file(templatefile)

        # Instantiate the cell from the template
        print("Loading cell %s from template file" % templatename)
        template = getattr(self.caller.hocObj, templatename)
     
        self.setCell(template(int(synapses)))
        print("exit2")

    def getLoadedTemplate(self, templatename):
        self.setCell(getattr(h, templatename)[0])

    def setCell(self,toWhat):
        self.cell = toWhat
        self.caller.cell = toWhat

    def get_main_ax2(self):
        main_ax = self.caller.hocObj.SectionList()
        at_terminal = 0

        current_secref = self.caller.hocObj.SectionRef(sec = self.caller.CellLoader.cell.axon[0])

        while not at_terminal:
            main_ax.append(current_secref.sec)
            if current_secref.nchild() > 0 :
                biggest_branch_diam = 0
                biggest_branch_ind = 0
                for i in range(current_secref.nchild()):
                    if current_secref.child[i](0).diam > biggest_branch_diam:
                        biggest_branch_diam = current_secref.child[i](0).diam
                        biggest_branch_ind = i
                current_secref = self.caller.hocObj.SectionRef(sec = current_secref.child[biggest_branch_ind])
            else:
                at_terminal = 1
            
                terminal_sec_str = self.caller.hocObj.secname(sec = current_secref.sec)

      
                for i in range(self.caller.interpCoordinates.numSect):
                    current_sec_str = self.caller.hocObj.secname(sec = self.caller.interpCoordinates.secrefs[i].sec)
                    if(current_sec_str == terminal_sec_str):
                        min_sec_ind = i
                        if self.debug:
                            print("Terminal section: {}, index: {}\n".format(terminal_sec_str,min_sec_ind))
        return main_ax

class myelinBiophysics():

    debug = True

    def __init__(self, caller):
        self.caller = caller
    def myelin_biophys(self):
        axon_bp = self.get_axon_biophys()

        for sec in self.caller.CellLoader.cell.somatic:
        
            sec.cm = 1
        for sec in self.caller.CellLoader.cell.apical:
          
            sec.cm = 2
        for sec in self.caller.CellLoader.cell.basal:
          
            sec.cm = 2

        for sec in self.caller.Morphology.axonal:
           
            sec.insert('pas')
            sec.e_pas = axon_bp.x[12]
            sec.Ra = 100
            sec.cm = 1
            sec.g_pas = axon_bp.x[13]


        for sec in self.caller.Morphology.Node_secList:
            
            # insert Ca_HVA
            sec.insert("SKv3_1")
            # insert SK_E2
            # insert CaDynamics_E2
            sec.insert("Nap_Et2")
            sec.insert("K_Pst")
            sec.insert("K_Tst")
            # insert Ca_LVAst
            sec.insert("NaTa_t")           
            # insert NaTa2_t
       
            for x in sec:

                if x > 0 and x < 1:
                    if (self.caller.hocObj.ismembrane("NaTa_t", sec = sec)): sec(x).gNaTa_tbar_NaTa_t = 2*axon_bp.x[0]

                    if (self.caller.hocObj.ismembrane("K_Tst", sec = sec)):  sec(x).gK_Tstbar_K_Tst = axon_bp.x[1]
                    if (self.caller.hocObj.ismembrane("CaDynamics_E2", sec = sec)):   sec(x).gamma_CaDynamics_E2 = axon_bp.x[2]
                    if (self.caller.hocObj.ismembrane("Nap_Et2", sec = sec)):  sec(x).gNap_Et2bar_Nap_Et2 = axon_bp.x[3]
                    if (self.caller.hocObj.ismembrane("SK_E2", sec = sec)):  sec(x).gSK_E2bar_SK_E2 = axon_bp.x[4]
                    if (self.caller.hocObj.ismembrane("Ca_HVA", sec = sec)):  sec(x).gCa_HVAbar_Ca_HVA = axon_bp.x[5]
                    if (self.caller.hocObj.ismembrane("K_Pst", sec = sec)):  sec(x).gK_Pstbar_K_Pst = axon_bp.x[6]
                    if (self.caller.hocObj.ismembrane("SKv3_1", sec = sec)):   sec(x).gSKv3_1bar_SKv3_1 = axon_bp.x[7]
                    if (self.caller.hocObj.ismembrane("CaDynamics_E2", sec = sec)):  sec(x).decay_CaDynamics_E2 = axon_bp.x[8]
                    if (self.caller.hocObj.ismembrane("Ca_LVAst", sec = sec)):  sec(x).gCa_LVAstbar_Ca_LVAst = axon_bp.x[9]
                    if (self.caller.hocObj.ismembrane("Im", sec = sec)):  sec(x).gImbar_Im = axon_bp.x[10]
                    if (self.caller.hocObj.ismembrane("Ca", sec = sec)): sec(x).gCabar_Ca = axon_bp.x[11]

            sec.ena = 50
            sec.ek = -85


    def get_axon_biophys(self):
        axon_bp = self.caller.hocObj.Vector(14)
        sec = self.caller.CellLoader.cell.axon[0]

        if (self.caller.hocObj.ismembrane("NaTa_t", sec = sec)): axon_bp.x[0] = sec.gNaTa_tbar_NaTa_t
        if (self.caller.hocObj.ismembrane("K_Tst", sec = sec)):  axon_bp.x[1] = sec.gK_Tstbar_K_Tst
        if (self.caller.hocObj.ismembrane("CaDynamics_E2", sec = sec)):  axon_bp.x[2] = sec.gamma_CaDynamics_E2
        if (self.caller.hocObj.ismembrane("Nap_Et2", sec = sec)):  axon_bp.x[3] = sec.gNap_Et2bar_Nap_Et2
        if (self.caller.hocObj.ismembrane("SK_E2", sec = sec)):  axon_bp.x[4] = sec.gSK_E2bar_SK_E2
        if (self.caller.hocObj.ismembrane("Ca_HVA", sec = sec)):  axon_bp.x[5] = sec.gCa_HVAbar_Ca_HVA
        if (self.caller.hocObj.ismembrane("K_Pst", sec = sec)):  axon_bp.x[6] = sec.gK_Pstbar_K_Pst
        if (self.caller.hocObj.ismembrane("SKv3_1", sec = sec)):  axon_bp.x[7] = sec.gSKv3_1bar_SKv3_1
        if (self.caller.hocObj.ismembrane("CaDynamics_E2", sec = sec)):  axon_bp.x[8] = sec.decay_CaDynamics_E2
        if (self.caller.hocObj.ismembrane("Ca_LVAst", sec = sec)):  axon_bp.x[9] = sec.gCa_LVAstbar_Ca_LVAst
        if (self.caller.hocObj.ismembrane("Im", sec = sec)):  axon_bp.x[10] = sec.gImbar_Im
        if (self.caller.hocObj.ismembrane("Ca", sec = sec)):  axon_bp.x[11] = sec.gCabar_Ca
        axon_bp.x[12] = sec.e_pas
        axon_bp.x[13] = sec.g_pas
        return axon_bp


class Morphology():

    debug = True

    def debug(self):
        if self.debug:
            pdb.set_trace()

    def __init__(self, caller):

        self.caller = caller

        self.INL_ratio = refpointer("INL_ratio", "100", hocObj = self.caller.hocObj)
        self.INL_ratio_term = refpointer("INL_ratio_term", "70", hocObj = self.caller.hocObj)
        self.nodeL = refpointer("nodeL", "1", hocObj = self.caller.hocObj)

        self.min_myelinL = refpointer("min_myelinL", "20", hocObj = self.caller.hocObj)
        self.min_myelinD = refpointer("min_myelinD", "0.2", hocObj = self.caller.hocObj)

        self.min_PMAS = refpointer("min_PMAS", "50", hocObj = self.caller.hocObj)

        self.myelinL_error = refpointer("myelinL_error", "0.1", hocObj = self.caller.hocObj)
        self.nodeL_error = refpointer("nodeL_error", "0.1", hocObj = self.caller.hocObj)
        self.max_myelin_order = refpointer("max_myelin_order", "0", hocObj = self.caller.hocObj)
        

        self.iseg_secList = objref("iseg_secList", hocObj = self.caller.hocObj)
        self.Node_secList = objref("Node_secList", hocObj = self.caller.hocObj)
        self.Myelin_secList = objref("Myelin_secList", hocObj = self.caller.hocObj)
        self.Unmyelin_secList = objref("Unmyelin_secList", hocObj = self.caller.hocObj)
        self.axonal = objref("axonal", hocObj = self.caller.hocObj)
        self.myelinCnts = objref("myelinCnts", hocObj = self.caller.hocObj)

        self.myelinBiophysics = myelinBiophysics(self.caller)

    def scale_diam2(self, factor, scale_seclist):

        diams = self.caller.hocObj.Vector()
        for sec in scale_seclist:
            diam_sec = self.getpt3d(5, sec=sec)
            diams.append(diam_sec)

        diams2 = diams.mul(factor)

        i = 0
        for sec in scale_seclist:
            for ii in range(int(self.caller.hocObj.n3d(sec=sec))):
                self.caller.hocObj.pt3dchange(ii, diams2.x[i])
                i = i + 1

    def scale_diam3(self, seclist):
        diams = self.caller.hocObj.Vector()
        g_ratios_sec = self.caller.hocObj.Vector()

        p1 = 0.6425 # polynomial curve fit of d vs. g_ratio from Micheva 2016 data
        p2 = -1.778
        p3 = 1.749
        p4 = 0.1518

        g_ratios = self.caller.hocObj.Vector()

        for sec in seclist:
            diams_sec = self.getpt3d(5, sec)
            diams.append(diams_sec)
           
            g_ratios_sec = diams_sec.c().pow(3).mul(p1)
            g_ratios_sec = g_ratios_sec.c().add(diams_sec.c().pow(2).mul(p2))
            g_ratios_sec = g_ratios_sec.c().add(diams_sec.c().mul(p3))
            g_ratios_sec = g_ratios_sec.c().add(p4)

            for nn in range(g_ratios_sec.size()):
                if g_ratios_sec[nn] > 0.8: g_ratios_sec[nn] = 0.8
                if g_ratios_sec[nn] < 0.4: g_ratios_sec[nn] = 0.4
            g_ratios.append(g_ratios_sec)
        ones = self.caller.hocObj.Vector(g_ratios.size())
        ones = ones.c().fill(1)

        diams2 = diams.c().mul(ones.c().div(g_ratios))
        i=0
        for sec in seclist:
            for ii in range(int(self.caller.hocObj.n3d(sec=sec))):
                self.caller.hocObj.pt3dchange(ii, diams2.x[i]) 
                i = i + 1



    def myelinate_axon(self, axon_secList):
        self.add_myelin(axon_secList)

        self.scale_diam3(axon_secList)
        
        self.geom_nseg(self.axonal, 40)

        self.myelinBiophysics.myelin_biophys()
        print("BioPhysics added")
        for sec in self.axonal:
            self.caller.CellLoader.cell.all.append(sec=sec)
            self.caller.CellLoader.cell.axonal.append(sec=sec)

    

        for sec in self.iseg_secList:
            self.caller.CellLoader.cell.axonal.append(sec=sec)

         

        for sec in self.axonal:
            self.caller.CellLoader.cell.axonal.append(sec=sec)
        print("sections added")
    

    def add_myelin(self, axon_secList):


        self.iseg_secList = self.caller.hocObj.SectionList()
        self.Myelin_secList = self.caller.hocObj.SectionList()
        self.Node_secList = self.caller.hocObj.SectionList()
        self.Unmyelin_secList = self.caller.hocObj.SectionList()
        self.axonal = self.caller.hocObj.SectionList()
        

        self.iseg_secList.append(self.caller.CellLoader.cell.axon[0])

        axon_secList.remove(self.iseg_secList)

        for sec in self.iseg_secList:
    
            if sec.L >= self.min_PMAS[0] + self.min_myelinL[0]:
                numMyelin = 1
                numNode = 1
                include_PMAS_myelin = 1
            else:
                numMyelin = 0
                numNode = 0
                include_PMAS_myelin = 0

        numAxonal = 0
        numUnmyelin = 0
        self.myelinCnts = self.caller.hocObj.Vector()

        max_order = self.get_max_order(axon_secList)

        for sec in axon_secList:
            # print(sec, sec(0).type_xtra, sec(1).type_xtra)
            numAxonal = numAxonal + 1
            if sec(1).type_xtra==2 or sec(1).type_xtra==5 :
                myelinL = sec(0).diam*self.INL_ratio_term[0]
            else:
                myelinL = sec(0).diam*self.INL_ratio[0]

            if sec(0).diam >= self.min_myelinD[0]:

                
                if self.check_in_secList(self.caller.CellLoader.main_ax_list, sec) \
                  or sec.order_xtra < max_order + 1 - self.max_myelin_order[0]:
                  
                    numMyelin_sec = int(sec.L/(myelinL + self.nodeL[0]))
                    # print("sec: {} L: {}, result: {}".format(sec, sec.L,  numMyelin_sec))
                    if numMyelin_sec == 0:
                       numMyelin_sec = int(sec.L/(self.min_myelinL[0] + self.nodeL[0]))
                else:
                    numMyelin_sec = 0   
            else:
                numMyelin_sec = 0
            
            numMyelin = numMyelin + numMyelin_sec
            numNode = numNode + numMyelin_sec
            if numMyelin_sec == 0:
                numUnmyelin = numUnmyelin + 1

            self.myelinCnts.append(numMyelin_sec)

        self.Myelin = create("Myelin", size = str(numMyelin),  hocObj = self.caller.hocObj)
        self.Node = create("Node", size = str(numNode), hocObj = self.caller.hocObj)

        if numUnmyelin >=1:
            self.Unmyelin = create("Unmyelin", size = str(numUnmyelin), hocObj = self.caller.hocObj) 

        for sec in self.Myelin:
            self.Myelin_secList.append(sec=sec)
            self.axonal.append(sec=sec)
        for sec in self.Node:
            self.Node_secList.append(sec=sec)
            self.axonal.append(sec=sec)
        for sec in self.Unmyelin:
            self.Unmyelin_secList.append(sec=sec)
            self.axonal.append(sec=sec)
 
        print("Myelinating axon: Replacing {} Axonal sections w/ {} Myelin, {} Node, {} Unmyelin sections\n".format(numAxonal,numMyelin,numNode,numUnmyelin))
        
        if include_PMAS_myelin == 1:
            print("Adding myelin before the 1st bifurcation")
            for sec in self.iseg_secList:

                children_SecList = self.getchildren(sec)
             
                for csec in children_SecList:
                    self.caller.hocObj.disconnect(sec=csec)
                self.Myelin[0].connect(sec(1), 0)

                secx = self.getpt3d(1,sec=sec)
                secy = self.getpt3d(2, sec=sec)
                secz = self.getpt3d(3, sec=sec)
                length = self.getpt3d(4, sec=sec)
                diamvec = self.getpt3d(5, sec=sec)
                last_pt3d_i = length.indwhere(">", self.min_PMAS[0]) -1 # is this correct conversion?

                while self.caller.hocObj.arc3d(self.caller.hocObj.n3d(sec=sec)-1, sec=sec) > self.min_PMAS[0]:
                    self.caller.hocObj.pt3dremove(self.caller.hocObj.n3d(sec=sec)-1, sec=sec)

                # get myelinL
                myelinL = length.x[-1] - self.min_PMAS[0] - self.nodeL[0]

            sec = self.Myelin[0]
            first_pt3d_i = last_pt3d_i
            last_pt3d_i = length.indwhere(">",self.min_PMAS[0] + myelinL)
            if last_pt3d_i > first_pt3d_i:
                last_pt3d_i = last_pt3d_i -1
            last_pt3d_i = self.add_new_points(first_pt3d_i, last_pt3d_i,secx,secy,secz,length,diamvec,myelinL, self.myelinL_error[0], 1, sec)
            
            self.Node[0].connect(self.Myelin[0](1), 0)

            sec = self.Node[0]
            first_pt3d_i = last_pt3d_i
            last_pt3d_i = length.size()-1
            if last_pt3d_i >  first_pt3d_i:
                last_pt3d_i = last_pt3d_i - 1
            else:
                last_pt3d_i = first_pt3d_i

            last_pt3d_i = self.add_new_points(first_pt3d_i,last_pt3d_i,secx,secy,secz,length,diamvec,self.nodeL[0],self.nodeL_error[0],1,sec)

            for csec in children_SecList:
                childe_SecRef = self.caller.hocObj.SectionRef(sec=csec)
                
                self.caller.hocObj.disconnect(sec=csec)
                childe_SecRef.sec.connect(self.Node[0](1), 0)

            mye_cnt = 1

        else:
            mye_cnt = 0
            print("No myelin before 1st bifurcation")

        sec_cnt = 0
        unmye_cnt = 0

        for sec in axon_secList:

            # print("________" + self.caller.hocObj.secname(sec=sec) + "__________")

            secx = self.getpt3d(1, sec=sec)
            secy = self.getpt3d(2, sec=sec)
            secz = self.getpt3d(3, sec=sec)
            length = self.getpt3d(4, sec=sec)
            diamvec = self.getpt3d(5, sec=sec)

            children_SecList = self.getchildren(sec)
            parent_SecList = self.getparent(sec)
            numMyelin_sec = int(self.myelinCnts.x[sec_cnt])

            myelinL_sec = 0
            self.caller.hocObj.delete_section(sec=sec)  # deleting sec?

            if numMyelin_sec >= 1:
                for parentsec in parent_SecList:
                    
                    self.Myelin[mye_cnt].connect(parentsec(1), 0)

                    secx.x[0] = self.caller.hocObj.x3d(self.caller.hocObj.n3d(sec=parentsec) - 1, sec=parentsec)
                    secy.x[0] = self.caller.hocObj.y3d(self.caller.hocObj.n3d(sec=parentsec) - 1, sec=parentsec)
                    secz.x[0] = self.caller.hocObj.z3d(self.caller.hocObj.n3d(sec=parentsec) - 1, sec=parentsec)
                    diamvec.x[0] = self.caller.hocObj.diam3d(self.caller.hocObj.n3d(sec=parentsec) - 1, sec=parentsec)
                    length = self.get_arc3d(secx, secy, secz)

                    # print(secx.size(), secy.size(), secz.size(), diamvec.size(), sec)

                    myelinL = (length.x[length.size()-1] - numMyelin_sec*self.nodeL[0])/numMyelin_sec
                # print("connect 5", self.Node[mye_cnt](0), self.Myelin[mye_cnt](1))
                self.Node[mye_cnt].connect(self.Myelin[mye_cnt](1), 0)
                
                # print("mL = {} mylenL = {}".format(self.Myelin[mye_cnt].L, myelinL_sec, ))

                first_pt3d_i = 0
                last_pt3d_i = length.indwhere(">", myelinL)
                if last_pt3d_i > first_pt3d_i:
                    last_pt3d_i = last_pt3d_i - 1
                # print("Add new points 1")
                last_pt3d_i = self.add_new_points(first_pt3d_i,last_pt3d_i,secx,secy,secz,length,diamvec,myelinL,self.myelinL_error[0],1, self.Myelin[mye_cnt])
                
                myelinL_sec = self.Myelin[mye_cnt].L

                # print("mL = {} mylenL = {}".format(self.Myelin[mye_cnt].L, myelinL_sec, ))

                if numMyelin_sec >= 2:
                    for nn in range(numMyelin_sec):

                        # print("L = {} mylenL = {}".format(self.Node[mye_cnt + nn].L, myelinL_sec))
                        # print("sec {}, myelinL {}, cnt: {}".format(h.secname(sec=self.Node[mye_cnt+nn]), myelinL_sec, mye_cnt))
                        first_pt3d_i = last_pt3d_i 
                       
                        last_pt3d_i = length.indwhere(">",myelinL_sec + self.nodeL[0])
                        # print(last_pt3d_i)
                        if last_pt3d_i > first_pt3d_i:
                            last_pt3d_i -= 1 
                        else: 
                            last_pt3d_i = first_pt3d_i
                        # print("Add new points 2")
                        last_pt3d_i = self.add_new_points(first_pt3d_i,last_pt3d_i,secx,secy,secz,length,diamvec,self.nodeL[0],self.nodeL_error[0],-1, self.Node[mye_cnt+nn])
                        
                        myelinL_sec = myelinL_sec + self.Node[mye_cnt + nn].L
                        if nn < numMyelin_sec - 1:
                            # print("connect 6")
                            self.Myelin[mye_cnt + nn + 1].connect(self.Node[mye_cnt + nn](1), 0)

                       

                        #  print("L = {} mylenL = {}".format(self.Node[mye_cnt + nn].L, myelinL_sec))


                        if nn < numMyelin_sec - 1:
                            # print("L = {} mylenL = {}".format(self.Myelin[mye_cnt + nn + 1].L, myelinL_sec))
                            first_pt3d_i = last_pt3d_i
                            last_pt3d_i = length.indwhere(">", myelinL_sec+myelinL)

                            if last_pt3d_i > first_pt3d_i:
                                last_pt3d_i = last_pt3d_i - 1
                            elif last_pt3d_i < 0: 
                                myelinL = length.x[length.size()-1] - length.x[first_pt3d_i] - self.nodeL[0]
                                last_pt3d_i = length.size() - 2
                            # print("Add new points 3")
                            last_pt3d_i =self.add_new_points(first_pt3d_i, last_pt3d_i, secx,secy, secz, length, diamvec, myelinL, self.myelinL_error[0], 1, self.Myelin[mye_cnt+nn+1])
                            
                            # print("connect 7")
                            myelinL_sec = myelinL_sec + self.Myelin[mye_cnt+nn+1].L
                            self.Node[mye_cnt+nn+1].connect(self.Myelin[mye_cnt+nn+1](1), 0)

                            # print("L = {} mylenL = {}".format(self.Myelin[mye_cnt + nn + 1].L, myelinL_sec))
                
                else:

                    first_pt3d_i = last_pt3d_i
                    last_pt3d_i = length.size() - 1
                    # print(self.caller.hocObj.secname(sec = self.Node[mye_cnt]), first_pt3d_i, last_pt3d_i)
                    # print("Add new points 4")
                    last_pt3d_i = self.add_new_points(first_pt3d_i,last_pt3d_i,secx,secy,secz,length,diamvec,self.nodeL[0],self.nodeL_error[0],1, self.Node[mye_cnt])


                for csec in children_SecList:

                    childe_SecRef = self.caller.hocObj.SectionRef(sec=csec)
                    # print("connect 8")
                    self.caller.hocObj.disconnect(sec=csec)
                    
                    childe_SecRef.sec.connect(self.Node[int(mye_cnt + numMyelin_sec - 1)](1), 0)
                
                mye_cnt += numMyelin_sec
           
            else:

                for psec in parent_SecList:
                    # print("connect 9", psec(1), self.Unmyelin[unmye_cnt](0))
                    # print("connecting: " + self.caller.hocObj.secname(sec=psec))

                    self.Unmyelin[unmye_cnt].connect(psec(1), 0)

                self.assign_pts(0,secx.size()-1, secx,secy,secz,diamvec, self.Unmyelin[unmye_cnt])

                for csec in children_SecList:
                    # print("disconnecting: " + self.caller.hocObj.secname(sec=csec))
                    childe_SecRef = self.caller.hocObj.SectionRef(sec=csec)
                    # print("connect 10")

                    self.caller.hocObj.disconnect(sec=csec)
                    # print("Connecting {} to {}".format(childe_SecRef.sec,self.Unmyelin[unmye_cnt](1) ))

                    childe_SecRef.sec.connect(self.Unmyelin[unmye_cnt](1), 0)

                unmye_cnt = unmye_cnt  + 1

            sec_cnt = sec_cnt + 1

    def getparent(self,sec):
        current_sec = self.caller.hocObj.SectionRef(sec=sec)
        parent = self.caller.hocObj.SectionList()

        parent.append(current_sec.parent)
        # print("Parent: " + self.caller.hocObj.secname(sec = current_sec.parent))
        return parent
 
    def getpt3d(self, dim, sec):
        nn = int(self.caller.hocObj.n3d(sec=sec))
        vec = self.caller.hocObj.Vector(nn)

        if (dim == 1):
            for ii in range(nn): 
                vec.x[ii] = self.caller.hocObj.x3d(ii, sec=sec)        
        
        elif(dim == 2):
            for ii in range(nn):
                vec.x[ii] = self.caller.hocObj.y3d(ii, sec=sec) 

        elif (dim == 3):
            for ii in range(nn):
                vec.x[ii] = self.caller.hocObj.z3d(ii, sec=sec)  

        elif (dim == 4):
            for ii in range(nn):
                 vec.x[ii] = self.caller.hocObj.arc3d(ii, sec=sec)  
            

        elif (dim == 5):
            for ii in range(nn):
                vec.x[ii] = self.caller.hocObj.diam3d(ii, sec=sec)
        

        return vec





    def add_new_points(self, first_pt3d_i, last_pt3d_i,secx,secy,secz,length,diamvec,secL, secerr, dir, sec):
        # print(length.x[last_pt3d_i] , length.x[first_pt3d_i] , secL , secerr)
        #print("add_new_points args: ", first_pt3d_i, last_pt3d_i,secx,secy,secz,length,diamvec,secL, secerr, dir, sec)
        # for (x,y,z,l, d) in zip(secx.x, secy.x, secz.x, length.x, diamvec.x):
        #   print("Vector contents: ",  x,y,z,l, d)


        if length.x[last_pt3d_i] - length.x[first_pt3d_i] >= secL + secerr:
            # print("dist1: {} - {} = {}, secL = {}, secerr = {}\n".format( last_pt3d_i,first_pt3d_i,length.x[last_pt3d_i] - length.x[first_pt3d_i],secL,secerr))
            last_pt3d_i = last_pt3d_i - 1
            self.assign_pts(first_pt3d_i,last_pt3d_i,secx,secy,secz, diamvec, sec)

            secL_add = secL - length.x[last_pt3d_i] - length.x[first_pt3d_i]

            interp_pt = self.add_interp_pt(last_pt3d_i, secx,secy,secz,diamvec,secL_add,dir, sec)

            """
            secx.x[last_pt3d_i] = interp_pt.x[0]
            secy.x[last_pt3d_i] = interp_pt.x[1]
            secz.x[last_pt3d_i] = interp_pt.x[2]
            length.x[last_pt3d_i] = length.x[last_pt3d_i-1] + \
                                    np.sqrt( (secx.x[last_pt3d_i] - interp_pt.x[0])**2 + \
                                    (secy.x[last_pt3d_i]-interp_pt.x[1])**2 + \
                                    (secz.x[last_pt3d_i]-interp_pt.x[2])**2 )  

            """

            secx.insrt(last_pt3d_i+1,interp_pt.x[0])
            secy.insrt(last_pt3d_i+1,interp_pt.x[1])
            secz.insrt(last_pt3d_i+1,interp_pt.x[2])
            length.insrt(last_pt3d_i+1,length.x[last_pt3d_i] + np.sqrt((secx.x[last_pt3d_i] - interp_pt.x[0])**2 + (secy.x[last_pt3d_i]-interp_pt.x[1])**2 + (secz.x[last_pt3d_i]-interp_pt.x[2])**2) ) # insert distnace of new point from first point

            diamvec.insrt(last_pt3d_i+1,diamvec.x[last_pt3d_i])
            # print(last_pt3d_i)
            return last_pt3d_i + 1

        elif length.x[last_pt3d_i] - length.x[first_pt3d_i] <= (secL - secerr):
            

            self.assign_pts(first_pt3d_i,last_pt3d_i,secx,secy,secz, diamvec, sec)
            # print("assign_pts sec Length:", sec.L)

            secL_add = secL - (length.x[last_pt3d_i] - length.x[first_pt3d_i])
        
            interp_pt = self.add_interp_pt(last_pt3d_i, secx, secy, secz, diamvec, secL_add, dir, sec)
            
            # print("interp_pt sec Length:", sec.L)

            secx.insrt(last_pt3d_i+1,interp_pt.x[0])
            secy.insrt(last_pt3d_i+1,interp_pt.x[1])
            secz.insrt(last_pt3d_i+1,interp_pt.x[2])        
        
            length.insrt(last_pt3d_i+1, length.x[last_pt3d_i] + np.sqrt( (secx.x[last_pt3d_i] - interp_pt.x[0])**2 + (secy.x[last_pt3d_i]-interp_pt.x[1])**2 + (secz.x[last_pt3d_i]-interp_pt.x[2])**2 )) 
            
            # print("sec Length:", sec.L)

            diamvec.insrt(last_pt3d_i+1,diamvec.x[last_pt3d_i])
            return last_pt3d_i+1 # position of new point in modified coordinate vectors  

        else:
            # print("dist3: {} - {} = {}, secL = {}, secerr = {}\n".format( last_pt3d_i,first_pt3d_i,length.x[last_pt3d_i] - length.x[first_pt3d_i],secL,secerr))
            # print(first_pt3d_i, last_pt3d_i, secx.size(), secy.size(), secz.size(), diamvec.size(), sec)
            self.assign_pts(first_pt3d_i, last_pt3d_i, secx, secy, secz, diamvec, sec) 
            # print("else - assign_pts sec Length:", sec.L)
            return last_pt3d_i

    # sets 2 comp per chunkSize in specified sectionList
    # geom_nseg(chunkSize,secList)
    def geom_nseg(self, seclist, chunksize=40):
        secIndex = 0
        for sec in seclist:
            sec.nseg = 1 + 2 * int(sec.L / chunksize)
            secIndex += 1
            #print(sec, sec.nseg, "set")

    def assign_pts(self, i1, i2, x, y, z, diamvec, sec):

        for i in range(i1, i2 + 1):
            self.caller.hocObj.pt3dadd(x.x[i], y.x[i], z.x[i], diamvec.x[i], sec=sec)
    
    def get_arc3d(self, x, y, z):

        length = self.caller.hocObj.Vector(x.size())
        length.x[0] = 0
        for i in range(1, x.size()):
            distance = np.sqrt((x.x[i] - x.x[i - 1]) ** 2
                                + (y.x[i] - y.x[i - 1]) ** 2
                                + (z.x[i] - z.x[i - 1]) ** 2)

            length.x[i] = length.x[i-1] + distance
        return length


    def add_interp_pt(self, i1, x, y, z, diams, len, dir,sec):
      
        inds = self.caller.hocObj.Vector()
        xtemp = x
        if dir == 1:
            inds.append(i1, i1 + 1)
            x = x.ind(inds)
            y = y.ind(inds)
            z = z.ind(inds)

            xu = x.x[1] - x.x[0]
            yu = y.x[1] - y.x[0]
            zu = z.x[1] - z.x[0]

            r = np.sqrt(xu ** 2 + yu ** 2 + zu ** 2)
            xn = x.x[0] + len * xu / r
            yn = y.x[0] + len * yu / r
            zn = z.x[0] + len * zu / r

        elif dir == -1:
            inds.append(i1 - 1, i1)
            x = x.ind(inds)
            y = y.ind(inds)
            z = z.ind(inds)

            xu = x.x[1] - x.x[0]
            yu = y.x[1] - y.x[0]
            zu = z.x[1] - z.x[0]

            r = np.sqrt(xu ** 2 + yu ** 2 + zu ** 2)
            xn = x.x[1] + len * xu / r
            yn = y.x[1] + len * yu / r
            zn = z.x[1] + len * zu / r 

        self.caller.hocObj.pt3dadd(xn, yn, zn, diams.x[i1], sec=sec)
        if dir == 1:
            dist = np.sqrt((xn - x.x[0]) ** 2 + (yn - y.x[0]) ** 2 + (zn - z.x[0]) ** 2)
        else:
            dist = np.sqrt( (xn - x.x[1])**2 + (yn - y.x[1])**2 + (zn - z.x[1])**2  )

        interp_pt = self.caller.hocObj.Vector()
        interp_pt.append(xn,yn,zn)


        return interp_pt

    # PROBABLY A BUG HERE I THINK, SHOULD WE REFENCE SEC?
    # yes - sectionlist.children function needs it
    def getchildren(self, sec):

        children = self.caller.hocObj.SectionList()
        children.children(sec = sec)
  
        return children


    def get_max_order(self, seclist):
        max_order = 0
        for sec in seclist:
            if self.caller.hocObj.ismembrane('xtra', sec=sec):
                if sec.order_xtra > max_order:
                    max_order = sec.order_xtra
            else:
                print("xtra not inserted in {}".format(self.caller.hocObj.secname(sec=sec)))
        return max_order

    def check_in_secList(self, seclist, currentsec):

        temp_seclist = self.caller.hocObj.SectionList()
        for sec in seclist:
            temp_seclist.append(sec)
        temp_seclist.append(sec=currentsec)
      
        return temp_seclist.unique() > 0




class MyelinatedCell():

    cellsCreated = 0

    def __init__(self, hocObj=None):

        if hocObj is None:
            self.hocObj = neuron.h
        else:
            self.hocObj = hocObj
        print("checking current secs")
        
        self.interpCoordinates = interpCoordinates(self)
        
        self.CellLoader = CellLoader(self)
        
        self.Morphology = Morphology(self)
        
        self.templatename   = None
        #createpanels()
        self.soma_point3 = objref("soma_point3", hocObj = self.hocObj)

    def color_plotmax(self, plot_mode = 1, save_fig = 0):
        self.hocObj.load_file(0, "anatscale.hoc")

    def loadcell(self, celllabel,myelinate_ax = True, loadedtemplate = None, synapses = True):
        

        if loadedtemplate is None:

            self.template = self.CellLoader.getTemplate('neurons/' + celllabel + '/')
        self.CellLoader.cell_chooser(celllabel, myelinate_ax=myelinate_ax, load_synapses=synapses, loadedtemplate  = loadedtemplate )
        self.cell = self.CellLoader.cell
        MyelinatedCell.cellsCreated += 1


        self.collect_hocobjects()
        print("objcollected")


    def collect_hocobjects(self):
        self.sections = self.hocObj.SectionList()
        for sec in self.Morphology.axonal:
            self.cell.all.append(sec=sec)
            self.cell.axonal.append(sec=sec)







def main(sim=False):
    neuron = 'L23_PC_cADpyr229_1'
    neuronfolder = 'neurons/' + neuron + '/'



    cwd = os.getcwd()
    os.chdir(neuronfolder)
    neuron.h.load_file("stdrun.hoc")
    neuron.h.load_file("import3d.hoc")


    neuron.h.load_file('stdlib.hoc')    #NEURON std. library
    neuron.h.load_file('import3d.hoc')  #import 3D morphology lib

    print('Loading constants')
    neuron.h.load_file('constants.hoc')

    neuron.h.load_file(1, "morphology.hoc")
    neuron.h.load_file(1, "biophysics.hoc")

    print('Loading constants')
    #neuron.h.load_file(1, 'template.hoc')
    neuron.h.load_file('synapses/synapses.hoc')

    celli = MyelinatedCell()

    celli.loadcell(neuron, myelinate_ax  = True, synapses=False)

    

    #init_simulation()
    #cell = create_cell()
    cell= celli.cell
    

    """
    stim_mode = 1 # 1 - ICMS, 2 - uniform E-field 
    xe = 200 # um electrode default position
    ye = -50 # um
    ze = 0  # um 
    sigma_e = 2.76e-7 # S/ubashm - conductivity in GM (Bungert 2016)
   
    # uniform E field stimulation
    theta = 180 # deg - polar angle
    phi = 0 # deg - azimuthal angle 

    amplitude = -1000
    
    Ex = amplitude * np.sin(theta) * np.cos(phi)
    Ey = amplitude * np.sin(theta) * np.sin(phi)
    Ez = amplitude * np.cos(theta)
    """
    
    for sec in neuron.h.allsec():
        if neuron.h.ismembrane("xtra", sec=sec):
            for x in sec:
                x.es_xtra = 0

    neuron.h._ref_stim_xtra[0] = 0

    """
    stim_amp = h.Vector(6)
    stim_time = h.Vector(6)

    stim_amp.fill(0)
    stim_amp.x[2] = 1
    stim_amp.x[3] = 1

    stim_time[1] = 1
    stim_time[2] = 1
    stim_time[3] = 4
    stim_time[4] = 4
    stim_time[5] = 5



    stim_amp.play(h._ref_stim_xtra, stim_time, 1)
    """
    if sim:

        stimuli = create_stimuli(cell, 1)

        recordings = {}

        recordings['time'] = neuron.h.Vector()
        recordings['soma(0.5)'] = neuron.h.Vector()

        recordings['time'].record(neuron.h._ref_t, 0.1)
        # recordings['soma(0.5)'].record(init.CellLoader.cell.soma[0](0.5)._ref_v, 0.1)
        recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_v, 0.1)

        neuron.h.tstop = 3000
        # create_stimuli(init.CellLoader.cell, 3)

        assert(sys.getrefcount(stimuli) > 1)

        print('Disabling variable timestep integration')
        neuron.h.cvode_active(0)
        neuron.h.finitialize(-65)
        neuron.h.fcurrent()
        print('Running for %f ms' % neuron.h.tstop)
        
        memireclist = neuron.h.List()
        for sec in neuron.h.allsec():
            for seg in sec:
                print(seg)
                memirec = neuron.h.Vector(int(neuron.h.tstop / neuron.h.dt+1))
                memirec.record(seg._ref_i_membrane, neuron.h.dt)
                memireclist.append(memirec)

        counter = 0

        neuron.h.t=0
        interval = 1000
        while neuron.h.t < neuron.h.tstop:
            
            neuron.h.fadvance()
            if counter == interval:
                print(neuron.h.t, recordings['soma(0.5)'].x[-1])
                counter = 0
            counter+=1

        imem = np.array(memireclist)

       


        #neuron.neuron.h.run()


        time = np.array(recordings['time'])
        soma_voltage = np.array(recordings['soma(0.5)'])


        soma_voltage_filename = 'soma_voltage_step.dat'

        np.savetxt(
                soma_voltage_filename,
                np.transpose(
                   np.vstack((
                        time,
                        soma_voltage))))

        print('Soma voltage for step saved to: %s'
              % soma_voltage_filename)

        import matplotlib.pyplot as plt
        plt.plot(recordings['time'], recordings['soma(0.5)'])
        plt.show()

    


def create_stimuli(cell, step_number):
    """Create the stimuli"""

    print('Attaching stimulus electrodes')

    stimuli = []
    step_amp = [0] * 3

    with open('current_amps.dat', 'r') as current_amps_file:
        first_line = current_amps_file.read().split('\n')[0].strip()
        hyp_amp, step_amp[0], step_amp[1], step_amp[2] = first_line.split(
            ' ')

    iclamp = neuron.h.IClamp(0.5, sec=cell.soma[0])
    iclamp.delay = 700
    iclamp.dur = 2000
    iclamp.amp = float(step_amp[step_number - 1])
    print('Setting up step current clamp: '
          'amp=%f nA, delay=%f ms, duration=%f ms' %
          (iclamp.amp, iclamp.delay, iclamp.dur))

    stimuli.append(iclamp)

    hyp_iclamp = neuron.h.IClamp(0.5, sec=cell.soma[0])
    hyp_iclamp.delay = 0
    hyp_iclamp.dur = 3000
    hyp_iclamp.amp = float(hyp_amp)
    print('Setting up hypamp current clamp: '
          'amp=%f nA, delay=%f ms, duration=%f ms' %
          (hyp_iclamp.amp, hyp_iclamp.delay, hyp_iclamp.dur))

    stimuli.append(hyp_iclamp)

    return stimuli

def create_cell(add_synapses=False):
    """Create the cell model"""
    # Load morphology

    neuron.h.load_file("morphology.hoc")
    # Load biophysics
    neuron.h.load_file("biophysics.hoc")
    # Load main cell template
    neuron.h.load_file("template.hoc")

    # Instantiate the cell from the template

    print("Loading cell cADpyr232_L5_TTPC2_8052133265")
    cell = neuron.h.cADpyr232_L5_TTPC2_8052133265(1 if add_synapses else 0)
    return cell


def init_simulation():
    """Initialise simulation environment"""

    neuron.h.load_file("stdrun.hoc")
    neuron.h.load_file("import3d.hoc")

    print('Loading constants')
    neuron.h.load_file('constants.hoc')




if __name__ == '__main__':
    if '-sim' in sys.argv:
        sim=True 
    else:
        sim=False
    main(sim=sim)