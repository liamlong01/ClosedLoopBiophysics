import os
import pdb

import neuron
h = neuron.h

import numpy as np

class strdef():

    debug = True

    def __init__(self, var):
        assert(type(var) == str)
        self.name = var

        # initialize objref in hoc
        try:
            self.create(var)

        except Exception as e:
            print(e)
            if self.debug:
                pdb.set_trace()

    def set(self, toWhat):
        h("{} = \"{}\"".format(self.name, toWhat))
        
        try:
            assert(getattr(h, self.name) == toWhat)
        except AssertionError as e:
            print(e)
            pdb.set_trace()

    def get(self):
        return getattr(h, self.name)

    def create(self, var):
        h("strdef %s" % var)


def create(var, size = None):
    assert(type(var) == str)
    assert(type(size) == str)
    if size is None:
        h("create {}".format(var))
    else:
        h("create {}[{}]".format(var, size))
    return getattr(h, var)

def refpointer(var, val):
    assert(type(var) == str)
    h("{} = {}" .format(var, val))
    return getattr(h, "_ref_{}".format(var))

def objref(var):
    assert(type(var) == str)
    h("objref {}".format(var))
    return getattr(h, var)



class interpCoordinates():
    debug = False
    def __init__(self, caller):
        self.caller = caller

        self.xx = objref("xx")
        self.yy = objref("yy")
        self.zz = objref("zz")
        self.length = objref("length")

        self.xint = objref("xint")
        self.yint = objref("yint")
        self.zint = objref("zint")
        self.range = objref("range")

        self.current_dir = strdef("current_dir")
        self.current_dir.set(h.getcwd())

        self.xList = h.List()
        self.yList = h.List()
        self.zList = h.List()

        self.numSect = refpointer("numSect", "1")
        self.numComp = refpointer("numSect", "1")


        self.secrefs = objref("secref")

    def getcoords(self):
        print("get coords.....")
        self.xList.remove_all()
        self.yList.remove_all()
        self.zList.remove_all()

        for sec in h.allsec():
            if h.ismembrane("xtra"):

                nn = int( h.n3d())
                self.xx = h.Vector(nn)
                self.yy = h.Vector(nn)
                self.zz = h.Vector(nn)
                self.length = h.Vector(nn)

                for ii in range(nn):
                    self.xx.x[ii] = h.x3d(ii)
                    self.yy.x[ii] = h.y3d(ii)
                    self.zz.x[ii] = h.z3d(ii)
                    self.length.x[ii] = h.arc3d(ii)


                self.length.div(self.length.x[nn-1])

                self.range = h.Vector(h.nseg + 2)
                self.range.indgen(1/h.nseg)
                self.range.sub(1/(2*h.nseg))
                self.range.x[0] = 0
                self.range.x[h.nseg + 1] = 1

                self.xint = h.Vector(h.nseg+2)
                self.yint = h.Vector(h.nseg+2)
                self.zint = h.Vector(h.nseg+2)

                self.xint.interpolate(self.range, self.length, self.xx)
                self.yint.interpolate(self.range, self.length, self.yy)
                self.zint.interpolate(self.range, self.length, self.zz)

                # original code x_xtra(xr) = .....
                # x_xtra is xtra mechanism of cas
                # the below should give pointer to that mechanism at same location
                for ii in range(1,h.nseg+1):
                    xr = self.range.x[ii]
                    sec(xr).x_xtra = self.xint.x[ii] # defined where????????????
                    sec(xr).y_xtra = self.xint.x[ii] # defined where????????????
                    sec(xr).z_xtra = self.xint.x[ii] # defined where????????????

                self.xint.remove(h.nseg+1)
                self.xint.remove(0)
                self.yint.remove(h.nseg+1)
                self.yint.remove(0)
                self.zint.remove(h.nseg+1)
                self.zint.remove(0)

                self.xList.append(self.xint)
                self.yList.append(self.yint)
                self.zList.append(self.zint)


        self.numSect = self.xList.count()
        print("numSect was %i\n" %  self.numSect) 



    

    def getSecRefs(self):
        self.getcoords()
        self.secrefs = h.List()
        secnum = 0
        self.numComp = 0  

        for sec in h.allsec():
            if h.ismembrane("xtra"):
                self.secrefs.append(h.SectionRef()) #SectionRef?????????
                secnum += 1
                self.numComp = self.numComp + h.nseg
        print("Created List of SectionRefs for %g sections\n" % secnum)

        self.assign_section_types()
        self.assign_order2()


        del secnum


    def traverse_tree(self, order, current_secref):

        current_secref.sec.order_xtra = order
       
 
        if current_secref.nchild() == 0:
            if self.debug:
                print("Reached terminal at {}, order is {}".format(h.secname(sec = current_secref.sec), 
                    current_secref.sec.order_xtra))

            #return to branch point
            not_branch = 1
            while not_branch:
               
                current_secref = h.SectionRef(sec = current_secref.parent)
                if current_secref.nchild() > 1:
                    not_branch = 0
                   
                else:
                    order = order - 1
               
        elif current_secref.nchild() == 1:
            if self.debug:
                print("1 Children in {}, order is {}".format(current_secref.nchild(), h.secname(sec = current_secref.sec), 
                    current_secref.sec.order_xtra))
            
          
            current_secref = h.SectionRef(sec = current_secref.child[0] )
       

            self.traverse_tree(order, current_secref)

        else:
            if self.debug:
                rint("{} Children in {}, order is {}".format(current_secref.nchild(),
                h.secname(sec=current_secref.sec),current_secref.sec.order_xtra) )


            order = order + 1
            for j in range(current_secref.nchild()): 
           
            

                self.traverse_tree(order, h.SectionRef(sec = current_secref.child[j]))

 

    def assign_order2(self):
        #  h.object_push(self.secrefs[0].sec)
        i = 1
        for n in range(self.secrefs[0].nchild()):

            print("children of soma:    %i"%n)
            oseci_secref = h.SectionRef(sec = self.secrefs[0].child[n])
            self.traverse_tree(i, oseci_secref)
     

    def assign_section_types(self):
        
        for i in range(self.numSect): 
            if self.secrefs[i].has_trueparent() == 0:
              
                self.secrefs[i].sec.type_xtra = 1
            else:
               
              

                parent_sec = h.SectionRef(sec = self.secrefs[i].parent().sec)
                parent_nchildren = parent_sec.nchild()

                Li = h.L

                if h.ismembrane("pas", sec=self.secrefs[i].sec):

                    Lambdai = self.Lambda(neuron.h.cas(), 1)*1000
                else:
                    Lambdai = h.L

                nn = h.nseg


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
                            self.secrefs[i].sec(1).type_xtra = 2

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


class cellChooser():

    debug = True

    def __init__(self, caller):
        self.caller = caller
        self.cell_id = refpointer("cell_id", "0")

        self.cell_names = objref("cell_names")
        self.cells = objref("cell")
        self.nil = objref("nil")

        self.main_ax_list = objref("main_ax_list")
        self.strobj = objref("strobj")

        self.secnames = objref("secnames")

        self.cell_names = h.List()

        self.current_dir = strdef("current_dir")
        self.current_dir.set("getcwd()")

        self.tstr2 = strdef("tstr2")
        self.cell_dir = strdef("cell_dir")
        self.createsim_file = strdef("createsim_file")

        for cell in ['L1_NGC_DA_bNAC219_1',
                     'L1_NGC_DA_bNAC219_2',
                     'L1_NGC_DA_bNAC219_3',
                     'L1_NGC_DA_bNAC219_4',
                     'L1_NGC_DA_bNAC219_5',
                     'L23_PC_cADpyr229_1',
                     'L23_PC_cADpyr229_2',
                     'L23_PC_cADpyr229_3',
                     'L23_PC_cADpyr229_4',
                     'L23_PC_cADpyr229_5',
                     'L4_LBC_cACint209_1',
                     'L4_LBC_cACint209_2',
                     'L4_LBC_cACint209_3',
                     'L4_LBC_cACint209_4',
                     'L4_LBC_cACint209_5',
                     'L5_TTPC2_cADpyr232_1',
                     'L5_TTPC2_cADpyr232_2',
                     'L5_TTPC2_cADpyr232_3',
                     'L5_TTPC2_cADpyr232_4',
                     'L5_TTPC2_cADpyr232_5',
                     'L6_TPC_L4_cADpyr231_1',
                     'L6_TPC_L4_cADpyr231_2',
                     'L6_TPC_L4_cADpyr231_3',
                     'L6_TPC_L4_cADpyr231_4',
                     'L6_TPC_L4_cADpyr231_5']:

            dir = 'cells/' + cell
            tstr1 = strdef(cell)
            tstr1.set(dir)
            self.cell_names.append(h.String(tstr1.get()))

    def cell_chooser(self, cell_id, myelinate_ax = True):
        self.cell_id = cell_id
        if cell_id > 0:
            h("forall delete_section()")
            self.get_cell()

            for sec in h.allsec():  
                # print("Inserting xtra into {}".format(h.secname(sec=sec)))

                sec.insert('xtra')
                sec.insert('extracellular')


            # Scaling
            scale_axon_diameter = 1.322 #TODO figure how to set this paramter for now default in Aberra code
            self.caller.editMorphology.scale_diam2(scale_axon_diameter, h.cell.axonal)


            self.caller.interpCoordinates.getSecRefs()

            meth1cells = h.Vector()
            meth1cells.append(6, 8, 20, 21, 22, 25)

            if not meth1cells.contains(cell_id):
                self.main_ax_list = self.get_main_ax2()
                print("got main_ax 2")
            else:
                self.main_ax_list = self.get_main_ax()


            if myelinate_ax:
                self.caller.editMorphology.myelinate_axon(h.cell.axonal)

              
                   

                numComp = 0

                for sec in h.allsec():
                    print(sec, ".nseg: ", sec.nseg)
                    sec.insert ("xtra")
                    sec.insert("extracellular")
                    numComp = sec.nseg + numComp
            

            print("Inserted xtra and extracellular in all {} compartments\n".format(numComp))
            self.setpointers()
            h.topology()

    def setpointers(self):
        self.caller.interpCoordinates.getSecRefs() 
        print("hello")
        for sec in h.allsec():


            if h.ismembrane("xtra", sec=sec) and h.ismembrane("extracellular", sec=sec):
                # print(sec)
                # should be a way to do this without resortin to this ugly str coding but fine for now
                neuron.h("for (x,0){ setpointer ex_xtra(x), e_extracellular(x) }")



    def get_cell(self):
        synapses_enabled = 0
        self.cell_dir.set(self.cell_names[self.cell_id - 1].s)
        os.chdir(self.cell_dir.get())
        h.load_file("createsimulation.hoc")
        h.create_cell(synapses_enabled)

    def get_main_ax2(self):
        main_ax = h.SectionList()
        at_terminal = 0

        current_secref = h.SectionRef(sec = h.cell.axon[0])

        while not at_terminal:
            main_ax.append(current_secref.sec)
            if current_secref.nchild() > 0 :
                biggest_branch_diam = 0
                biggest_branch_ind = 0
                for i in range(current_secref.nchild()):
                    if current_secref.child[i](0).diam > biggest_branch_diam:
                        biggest_branch_diam = current_secref.child[i](0).diam
                        biggest_branch_ind = i
                current_secref = h.SectionRef(sec = current_secref.child[biggest_branch_ind])
            else:
                at_terminal = 1
                terminal_sec_str = h.String()
                terminal_sec_str.s = h.secname(sec = current_secref.sec)

                current_sec_str = h.String()
                for i in range(self.caller.interpCoordinates.numSect):
                    current_sec_str.s = h.secname(sec = self.caller.interpCoordinates.secrefs[i].sec)
                    if(current_sec_str.s == terminal_sec_str.s):
                        min_sec_ind = i
                        if self.debug:
                            print("Terminal section: {}, index: {}\n".format(terminal_sec_str.s,min_sec_ind))
        
        return main_ax

class myelinBiophysics():

    debug = True

    def __init__(self, caller):
        self.caller = caller
    def myelin_biophys(self):
        axon_bp = self.get_axon_biophys()

        for sec in h.cell.somatic:
            sec.cm = 1
        for sec in h.cell.apical:
            sec.cm = 2
        for sec in h.cell.basal:
            sec.cm = 2

        for sec in h.cell.axonal:
            sec.insert('pas')
            sec.e_pas = axon_bp.x[12]
            sec.Ra = 100
            sec.cm = 1
            sec.g_pas = axon_bp.x[13]


        for sec in self.caller.editMorphology.Node_secList:
            
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
                    if (h.ismembrane("NaTa_t", sec = sec)): sec(x).gNaTa_tbar_NaTa_t = 2*axon_bp.x[0]

                    if (h.ismembrane("K_Tst", sec = sec)):  sec(x).gK_Tstbar_K_Tst = axon_bp.x[1]
                    if (h.ismembrane("CaDynamics_E2", sec = sec)):   sec(x).gamma_CaDynamics_E2 = axon_bp.x[2]
                    if (h.ismembrane("Nap_Et2", sec = sec)):  sec(x).gNap_Et2bar_Nap_Et2 = axon_bp.x[3]
                    if (h.ismembrane("SK_E2", sec = sec)):  sec(x).gSK_E2bar_SK_E2 = axon_bp.x[4]
                    if (h.ismembrane("Ca_HVA", sec = sec)):  sec(x).gCa_HVAbar_Ca_HVA = axon_bp.x[5]
                    if (h.ismembrane("K_Pst", sec = sec)):  sec(x).gK_Pstbar_K_Pst = axon_bp.x[6]
                    if (h.ismembrane("SKv3_1", sec = sec)):   sec(x).gSKv3_1bar_SKv3_1 = axon_bp.x[7]
                    if (h.ismembrane("CaDynamics_E2", sec = sec)):  sec(x).decay_CaDynamics_E2 = axon_bp.x[8]
                    if (h.ismembrane("Ca_LVAst", sec = sec)):  sec(x).gCa_LVAstbar_Ca_LVAst = axon_bp.x[9]
                    if (h.ismembrane("Im", sec = sec)):  sec(x).gImbar_Im = axon_bp.x[10]
                    if (h.ismembrane("Ca", sec = sec)): sec(x).gCabar_Ca = axon_bp.x[11]

            sec.ena = 50
            sec.ek = -85


    def get_axon_biophys(self):
        axon_bp = h.Vector(14)
        sec = h.cell.axon[0]

        if (h.ismembrane("NaTa_t", sec = sec)): axon_bp.x[0] = sec.gNaTa_tbar_NaTa_t
        if (h.ismembrane("K_Tst", sec = sec)):  axon_bp.x[1] = sec.gK_Tstbar_K_Tst
        if (h.ismembrane("CaDynamics_E2", sec = sec)):  axon_bp.x[2] = sec.gamma_CaDynamics_E2
        if (h.ismembrane("Nap_Et2", sec = sec)):  axon_bp.x[3] = sec.gNap_Et2bar_Nap_Et2
        if (h.ismembrane("SK_E2", sec = sec)):  axon_bp.x[4] = sec.gSK_E2bar_SK_E2
        if (h.ismembrane("Ca_HVA", sec = sec)):  axon_bp.x[5] = sec.gCa_HVAbar_Ca_HVA
        if (h.ismembrane("K_Pst", sec = sec)):  axon_bp.x[6] = sec.gK_Pstbar_K_Pst
        if (h.ismembrane("SKv3_1", sec = sec)):  axon_bp.x[7] = sec.gSKv3_1bar_SKv3_1
        if (h.ismembrane("CaDynamics_E2", sec = sec)):  axon_bp.x[8] = sec.decay_CaDynamics_E2
        if (h.ismembrane("Ca_LVAst", sec = sec)):  axon_bp.x[9] = sec.gCa_LVAstbar_Ca_LVAst
        if (h.ismembrane("Im", sec = sec)):  axon_bp.x[10] = sec.gImbar_Im
        if (h.ismembrane("Ca", sec = sec)):  axon_bp.x[11] = sec.gCabar_Ca
        axon_bp.x[12] = sec.e_pas
        axon_bp.x[13] = sec.g_pas
        return axon_bp


class editMorphology():

    debug = True

    def debug(self):
        if self.debug:
            pdb.set_trace()

    def __init__(self, caller):

        self.caller = caller

        self.INL_ratio = refpointer("INL_ratio", "100")
        self.INL_ratio_term = refpointer("INL_ratio_term", "70")
        self.nodeL = refpointer("nodeL", "1")

        self.min_myelinL = refpointer("min_myelinL", "20")
        self.min_myelinD = refpointer("min_myelinD", "0.2")

        self.min_PMAS = refpointer("min_PMAS", "50")

        self.myelinL_error = refpointer("myelinL_error", "0.1")
        self.nodeL_error = refpointer("nodeL_error", "0.1")
        self.max_myelin_order = refpointer("max_myelin_order", "0")
        self.Myelin = create("Myelin", size="2")
        self.Node = create("Node", size="2")
        self.Unmyelin = create("Unmyelin", size="2")

        self.iseg_secList = objref("iseg_secList")
        self.Node_secList = objref("Node_secList")
        self.Myelin_secList = objref("Myelin_secList")
        self.Unmyelin_secList = objref("Unmyelin_secList")
        self.axonal = objref("axonal")
        self.myelinCnts = objref("myelinCnts")

        self.myelinBiophysics = myelinBiophysics(self.caller)

    def scale_diam2(self, factor, scale_seclist):

        diams = h.Vector()
        for sec in scale_seclist:
            diam_sec = self.getpt3d(5, sec=sec)
            diams.append(diam_sec)

        diams2 = diams.mul(factor)

        i = 0
        for sec in scale_seclist:
            for ii in range(int(h.n3d(sec=sec))):
                h.pt3dchange(ii, diams2.x[i])
                i = i + 1

    def scale_diam3(self, seclist):
        diams = h.Vector()
        g_ratios_sec = h.Vector()

        p1 = 0.6425 # polynomial curve fit of d vs. g_ratio from Micheva 2016 data
        p2 = -1.778
        p3 = 1.749
        p4 = 0.1518

        g_ratios = h.Vector()

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
        ones = h.Vector(g_ratios.size())
        ones = ones.c().fill(1)

        diams2 = diams.c().mul(ones.c().div(g_ratios))
        i=0
        for sec in seclist:
            for ii in range(int(h.n3d(sec=sec))):
                h.pt3dchange(ii, diams2.x[i])
                i = i + 1



    def myelinate_axon(self, axon_secList):
        self.add_myelin(axon_secList)

        self.scale_diam3(axon_secList)

        self.geom_nseg(self.axonal, 40)

        self.myelinBiophysics.myelin_biophys()
        print("BioPhysics added")
        for sec in self.axonal:
            h.cell.all.append(sec=sec)
            h.cell.axonal.append(sec=sec)

    

        for sec in self.iseg_secList:
            h.cell.axonal.append(sec=sec)

         

        for sec in self.axonal:
            h.cell.axonal.append(sec=sec)
        print("sections added")
    

    def add_myelin(self, axon_secList):


        self.iseg_secList = h.SectionList()
        self.Myelin_secList = h.SectionList()
        self.Node_secList = h.SectionList()
        self.Unmyelin_secList = h.SectionList()
        self.axonal = h.SectionList()
        

        self.iseg_secList.append(h.cell.axon[0])

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
        self.myelinCnts = h.Vector()

        max_order = self.get_max_order(axon_secList)

        for sec in axon_secList:
            # print(sec, sec(0).type_xtra, sec(1).type_xtra)
            numAxonal = numAxonal + 1
            if sec(1).type_xtra==2 or sec(1).type_xtra==5 :
                myelinL = sec(0).diam*self.INL_ratio_term[0]
            else:
                myelinL = sec(0).diam*self.INL_ratio[0]

            if sec(0).diam >= self.min_myelinD[0]:

                
                if self.check_in_secList(self.caller.cellChooser.main_ax_list, sec) \
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

        self.Myelin = create("Myelin", size = str(numMyelin))
        self.Node = create("Node", size = str(numNode))

        if numUnmyelin >=1:
            self.Unmyelin = create("Unmyelin", size = str(numUnmyelin)) 

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
                print("connect 1")
                for csec in children_SecList:
                    h.disconnect(sec=csec)
                sec.connect(self.Myelin[0](0), 1)

                secx = self.getpt3d(1,sec=sec)
                secy = self.getpt3d(2, sec=sec)
                secz = self.getpt3d(3, sec=sec)
                length = self.getpt3d(4, sec=sec)
                diamvec = self.getpt3d(5, sec=sec)
                last_pt3d_i = length.indwhere(">", self.min_PMAS[0]) -1 # is this correct conversion?

                while arcrd(n3d()-1) > self.min_PMAS[0]:
                    pt3dremove(n3d()-1)

                # get myelinL
                myelinL = length.x[-1] - self.min_PMAS[0] - self.nodeL[0]

            sec = self.Myelin[0]
            first_pt3d_i = laexitst_pt3d_i
            last_pt3d_i = length.indwhere(">",self.min_PMAS[0] + myelinL)
            if last_pt3d_i > first_pt3d_i:
                last_pt3d_i = last_pt3d_i -1
            last_pt3d_i = self.add_new_points(first_pt3d_i, last_pt3d_i,secx,secy,secz,length,diamvec,myelinL, self.myelinL_error, 1, sec)
            print("connect 2")
            self.Node[0].connect(Myelin[0](1), 0)

            sec = Node[0]
            first_pt3d_i = last_pt3d_i
            last_pt3d_i = length.size()-1
            if last_pt3d_i >  first_pt3d_i:
                last_pt3d_i = last_pt3d_i - 1
            else:
                last_pt3d_i = first_pt3d_i

            last_pt3d_i = self.add_new_points(first_pt3d_i,last_pt3d_i,secx,secy,secz,length,diamvec,self.nodeL[0],self.nodeL_error[0],1,sec)

            for csec in children_SecList:
                childe_SecRef = h.SectionRef(sec=csec)
                print("connect 3")
                h.disconnect(sec=csec)
                childe_SecRef.sec.connect(self.Node[0](1), 0)

            mye_cnt = 1

        else:
            mye_cnt = 0
            print("No myelin before 1st bifurcation")

        sec_cnt = 0
        unmye_cnt = 0

        for sec in axon_secList:

            print("________" + h.secname(sec=sec) + "__________")

            secx = self.getpt3d(1, sec=sec)
            secy = self.getpt3d(2, sec=sec)
            secz = self.getpt3d(3, sec=sec)
            length = self.getpt3d(4, sec=sec)
            diamvec = self.getpt3d(5, sec=sec)

            children_SecList = self.getchildren(sec)
            parent_SecList = self.getparent(sec)
            numMyelin_sec = int(self.myelinCnts.x[sec_cnt])

            myelinL_sec = 0
            h.delete_section(sec=sec)  # deleting sec?

            if numMyelin_sec >= 1:
                for parentsec in parent_SecList:
                    print("connect 4")
                    self.Myelin[mye_cnt].connect(parentsec(1), 0)

                    secx.x[0] = h.x3d(h.n3d(sec=parentsec) - 1, sec=parentsec)
                    secy.x[0] = h.y3d(h.n3d(sec=parentsec) - 1, sec=parentsec)
                    secz.x[0] = h.z3d(h.n3d(sec=parentsec) - 1, sec=parentsec)
                    diamvec.x[0] = h.diam3d(h.n3d(sec=parentsec) - 1, sec=parentsec)
                    length = self.get_arc3d(secx, secy, secz)

                    print(secx.size(), secy.size(), secz.size(), diamvec.size(), sec)

                    myelinL = (length.x[length.size()-1] - numMyelin_sec*self.nodeL[0])/numMyelin_sec
                print("connect 5")
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

                        print("L = {} mylenL = {}".format(self.Node[mye_cnt + nn].L, myelinL_sec))
                        # print("sec {}, myelinL {}, cnt: {}".format(h.secname(sec=self.Node[mye_cnt+nn]), myelinL_sec, mye_cnt))
                        first_pt3d_i = last_pt3d_i 
                       
                        last_pt3d_i = length.indwhere(">",myelinL_sec + self.nodeL[0])
                        # print(last_pt3d_i)
                        if last_pt3d_i > first_pt3d_i:
                            last_pt3d_i -= 1 
                        else: 
                            last_pt3d_i = first_pt3d_i
                        print("Add new points 2")
                        last_pt3d_i = self.add_new_points(first_pt3d_i,last_pt3d_i,secx,secy,secz,length,diamvec,self.nodeL[0],self.nodeL_error[0],-1, self.Node[mye_cnt+nn])
                        
                        myelinL_sec = myelinL_sec + self.Node[mye_cnt + nn].L
                        if nn < numMyelin_sec - 1:
                            print("connect 6")
                            self.Node[mye_cnt + nn].connect(self.Myelin[mye_cnt + nn + 1](0), 1)

                       

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
                            print("Add new points 3")
                            last_pt3d_i =self.add_new_points(first_pt3d_i, last_pt3d_i, secx,secy, secz, length, diamvec, myelinL, self.myelinL_error[0], 1, self.Myelin[mye_cnt+nn+1])
                            
                            print("connect 7")
                            myelinL_sec = myelinL_sec + self.Myelin[mye_cnt+nn+1].L
                            self.Myelin[mye_cnt+nn+1].connect(self.Node[mye_cnt+nn+1](0), 1)

                            # print("L = {} mylenL = {}".format(self.Myelin[mye_cnt + nn + 1].L, myelinL_sec))
                
                else:

                    first_pt3d_i = last_pt3d_i
                    last_pt3d_i = length.size() - 1
                    # print(h.secname(sec = self.Node[mye_cnt]), first_pt3d_i, last_pt3d_i)
                    print("Add new points 4")
                    last_pt3d_i = self.add_new_points(first_pt3d_i,last_pt3d_i,secx,secy,secz,length,diamvec,self.nodeL[0],self.nodeL_error[0],1, self.Node[mye_cnt])


                for csec in children_SecList:

                    childe_SecRef = h.SectionRef(sec=csec)
                    print("connect 8")
                    h.disconnect(sec=csec)
                    
                    childe_SecRef.sec.connect(self.Node[int(mye_cnt + numMyelin_sec - 1)](1), 0)
                
                mye_cnt += numMyelin_sec
           
            else:

                for psec in parent_SecList:
                    print("connect 9", psec(1), self.Unmyelin[unmye_cnt](0))
                    # print("connecting: " + h.secname(sec=psec))

                    self.Unmyelin[unmye_cnt].connect(psec(1), 0)

                self.assign_pts(0,secx.size()-1, secx,secy,secz,diamvec, self.Unmyelin[unmye_cnt])

                for csec in children_SecList:
                    # print("disconnecting: " + h.secname(sec=csec))
                    childe_SecRef = h.SectionRef(sec=csec)
                    print("connect 10")

                    h.disconnect(sec=csec)
                    # print("Connecting {} to {}".format(childe_SecRef.sec,self.Unmyelin[unmye_cnt](1) ))

                    childe_SecRef.sec.connect(self.Unmyelin[unmye_cnt](1), 0)

                unmye_cnt = unmye_cnt  + 1

            sec_cnt = sec_cnt + 1

    def getparent(self,sec):
        current_sec = h.SectionRef(sec=sec)
        parent = h.SectionList()

        parent.append(current_sec.parent)
        print("Parent: " + h.secname(sec = current_sec.parent))
        return parent
 
    def getpt3d(self, dim, sec):
        nn = int(h.n3d(sec=sec))
        vec = h.Vector(nn)

        if (dim == 1):
            for ii in range(nn): 
                vec.x[ii] = h.x3d(ii, sec=sec)        
        
        elif(dim == 2):
            for ii in range(nn):
                vec.x[ii] = h.y3d(ii, sec=sec) 

        elif (dim == 3):
            for ii in range(nn):
                vec.x[ii] = h.z3d(ii, sec=sec)  

        elif (dim == 4):
            for ii in range(nn):
                 vec.x[ii] = h.arc3d(ii, sec=sec)  
            

        elif (dim == 5):
            for ii in range(nn):
                vec.x[ii] = h.diam3d(ii, sec=sec)
        

        return vec





    def add_new_points(self, first_pt3d_i, last_pt3d_i,secx,secy,secz,length,diamvec,secL, secerr, dir, sec):
        # print(length.x[last_pt3d_i] , length.x[first_pt3d_i] , secL , secerr)
        #print("add_new_points args: ", first_pt3d_i, last_pt3d_i,secx,secy,secz,length,diamvec,secL, secerr, dir, sec)
        # for (x,y,z,l, d) in zip(secx.x, secy.x, secz.x, length.x, diamvec.x):
        #   print("Vector contents: ",  x,y,z,l, d)


        if length.x[last_pt3d_i] - length.x[first_pt3d_i] >= secL + secerr:
            print("dist1: {} - {} = {}, secL = {}, secerr = {}\n".format( last_pt3d_i,first_pt3d_i,length.x[last_pt3d_i] - length.x[first_pt3d_i],secL,secerr))
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
            print(last_pt3d_i)
            return last_pt3d_i + 1

        elif length.x[last_pt3d_i] - length.x[first_pt3d_i] <= (secL - secerr):
            

            self.assign_pts(first_pt3d_i,last_pt3d_i,secx,secy,secz, diamvec, sec)
            print("assign_pts sec Length:", sec.L)

            secL_add = secL - (length.x[last_pt3d_i] - length.x[first_pt3d_i])
        
            interp_pt = self.add_interp_pt(last_pt3d_i, secx, secy, secz, diamvec, secL_add, dir, sec)
            
            print("interp_pt sec Length:", sec.L)

            secx.insrt(last_pt3d_i+1,interp_pt.x[0])
            secy.insrt(last_pt3d_i+1,interp_pt.x[1])
            secz.insrt(last_pt3d_i+1,interp_pt.x[2])        
        
            length.insrt(last_pt3d_i+1, length.x[last_pt3d_i] + np.sqrt( (secx.x[last_pt3d_i] - interp_pt.x[0])**2 + (secy.x[last_pt3d_i]-interp_pt.x[1])**2 + (secz.x[last_pt3d_i]-interp_pt.x[2])**2 )) 
            
            print("sec Length:", sec.L)

            diamvec.insrt(last_pt3d_i+1,diamvec.x[last_pt3d_i])
            return last_pt3d_i+1 # position of new point in modified coordinate vectors  

        else:
            # print("dist3: {} - {} = {}, secL = {}, secerr = {}\n".format( last_pt3d_i,first_pt3d_i,length.x[last_pt3d_i] - length.x[first_pt3d_i],secL,secerr))
            # print(first_pt3d_i, last_pt3d_i, secx.size(), secy.size(), secz.size(), diamvec.size(), sec)
            self.assign_pts(first_pt3d_i, last_pt3d_i, secx, secy, secz, diamvec, sec) 
            print("else - assign_pts sec Length:", sec.L)
            return last_pt3d_i

    # sets 2 comp per chunkSize in specified sectionList
    # geom_nseg(chunkSize,secList)
    def geom_nseg(self, seclist, chunksize=40):
        secIndex = 0
        for sec in seclist:
            sec.nseg = 1 + 2 * int(sec.L / chunksize)
            secIndex += 1
            print(sec, sec.nseg, "set")

    def assign_pts(self, i1, i2, x, y, z, diamvec, sec):

        for i in range(i1, i2 + 1):
            h.pt3dadd(x.x[i], y.x[i], z.x[i], diamvec.x[i], sec=sec)
    
    def get_arc3d(self, x, y, z):

        length = h.Vector(x.size())
        length.x[0] = 0
        for i in range(1, x.size()):
            distance = np.sqrt((x.x[i] - x.x[i - 1]) ** 2
                                + (y.x[i] - y.x[i - 1]) ** 2
                                + (z.x[i] - z.x[i - 1]) ** 2)

            length.x[i] = length.x[i-1] + distance
        return length


    def add_interp_pt(self, i1, x, y, z, diams, len, dir,sec):
      
        inds = h.Vector()
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
 


        h.pt3dadd(xn, yn, zn, diams.x[i1], sec=sec)
        if dir == 1:
            dist = np.sqrt((xn - x.x[0]) ** 2 + (yn - y.x[0]) ** 2 + (zn - z.x[0]) ** 2)
        else:
            dist = np.sqrt( (xn - x.x[1])**2 + (yn - y.x[1])**2 + (zn - z.x[1])**2  )

        interp_pt = h.Vector()
        interp_pt.append(xn,yn,zn)


        return interp_pt

    # PROBABLY A BUG HERE I THINK, SHOULD WE REFENCE SEC?
    # yes - sectionlist.children function needs it
    def getchildren(self, sec):

        children = h.SectionList()
        children.children(sec = sec)
  
        return children


    def get_max_order(self, seclist):
        max_order = 0
        for sec in seclist:
            if h.ismembrane('xtra', sec=sec):
                if sec.order_xtra > max_order:
                    max_order = sec.order_xtra
            else:
                print("xtra not inserted in {}".format(h.secname(sec=sec)))
        return max_order

    def check_in_secList(self, seclist, currentsec):

        temp_seclist = h.SectionList()
        for sec in seclist:
            temp_seclist.append(sec)
        temp_seclist.append(sec=currentsec)
      
        return temp_seclist.unique() > 0




class initialize():

    def __init__(self):
        
        h.load_file("nrngui.hoc")
        self.interpCoordinates = interpCoordinates(self)
        h.load_file("setPointers.hoc")
        h.load_file("calcVe.hoc")
        h.load_file("stimWaveform.hoc")
        self.cellChooser = cellChooser(self)
        h.load_file("setParams.hoc")
        self.editMorphology = editMorphology(self)
        h.load_file("rig.ses")

        #createpanels()

        self.soma_point3 = objref("soma_point3")

    def color_plotmax(self, plot_mode = 1, save_fig = 0):
        h.load_file(0, "anatscale.hoc")






if __name__ == '__main__':
    init = initialize()

    init.cellChooser.cell_chooser(1)
