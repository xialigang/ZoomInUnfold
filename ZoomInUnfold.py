from ROOT import TFile, TTree, TH1F, TH2F
import math
import numpy as np
from array import array


class ZoomInUnfold(object):
    def __init__(self, sig_filepath, sig_treename, data_filepath, data_treename, mtrue, mrec, xmin, xmax, nbins, nsplit=0, res=1, dotoymc=0, ntoys=-1):

        self.sig_filepath = sig_filepath # file path to signal mc
        self.sig_treename = sig_treename # tree name in signal mc

        self.data_filepath = data_filepath # file path to datanal data
        self.data_treename = data_treename # tree name in datanal data


        self.mtrue = mtrue # true variable name in signal mc root file
        self.mrec = mrec   # rec. variable name in data root file
        self.xmin = xmin # min value of the varialbe
        self.xmax = xmax # max value
        self.nbins = nbins # initial number of bins
        self.nsplit = nsplit
        self.res = res

        self.dotoymc = dotoymc
        self.ntoys = ntoys

        self.htrue = None
        self.hrec = None
        self.hout = None
        self.hVx = None
        self.hRx= None

        self.Ltrue = []
        self.Lrec = []
        self.Lout = []
        self.Vx = None
        self.Rx = None
        self.Lisplit = [] # a list to record the splitting positions, toy mc needs this

        self.debug = 0
        self.initbinning =[]
        self.binning =[]
        dx = (xmax-xmin)/nbins
        for i in range(nbins+1): # initialize the binning
            self.initbinning.append(xmin+i*dx)
    def get_binning(self,binning, Isplit):
        oldbinning = binning
        newbinning = []
        for i in range(len(oldbinning)):
            if i == Isplit-1:
                newbinning.append(oldbinning[i])
                newbinning.append(0.5*(float(oldbinning[i])+float(oldbinning[i+1])))
            else:
                newbinning.append(oldbinning[i])
        return newbinning
    def get_b(self, tree, binning, dotoymc=0):
        nbins = len(binning)
        array_binning = np.array(binning)
        #hrec = TH1F('hrec', '', nbins-1, array_binning)
        hrec = TH1F('', '', nbins-1, array_binning)
        Lrec = []
        for irec in range(nbins-1):
            Nrec = 0.
            if dotoymc:
                Nrec = tree.GetEntries('%s>%.2f && %s<%.2f' %(self.mrec, binning[irec], self.mrec, binning[irec+1]))
            else:
                Nrec = tree.GetEntries('%s>%.2f && %s<%.2f && %s>%.2f && %s<%.2f' %(self.mrec, binning[irec], self.mrec, binning[irec+1], self.mtrue, self.xmin, self.mtrue, self.xmax))
            hrec.SetBinContent(irec+1, Nrec)
            hrec.SetBinError(irec+1, math.sqrt(Nrec))
            Lrec.append([Nrec])
        array_Lrec = np.array(Lrec)
        if self.debug:
            print('We are in get_b(...)')
            print('binning =',binning)
            print('hrec =', hrec, hrec.Integral())
        return hrec, array_Lrec
    def get_htrue(self, tree, binning):
        nbins = len(binning)
        array_binning = np.array(binning)
        htrue = TH1F('', '', nbins-1, array_binning)
        Ltrue = []
        for itrue in range(nbins-1):
            Ntrue = tree.GetEntries('%s>%.2f && %s<%.2f && %s>%.2f && %s<%.2f' %(self.mtrue, binning[itrue], self.mtrue, binning[itrue+1], self.mrec, self.xmin, self.mrec, self.xmax))
            htrue.SetBinContent(itrue+1, Ntrue)
            htrue.SetBinError(itrue+1, math.sqrt(Ntrue))
            Ltrue.append(Ntrue)
        array_Ltrue = np.array(Ltrue)
        return htrue
    def get_response_matrix(self, tree, binning):
        nbins = len(binning)
        A = []
        for itrue in range(nbins-1):
            Ntrue = tree.GetEntries('%s>%.2f && %s<%.2f && %s>%.2f && %s<%.2f' %(self.mtrue, binning[itrue], self.mtrue, binning[itrue+1], self.mrec, self.xmin, self.mrec, self.xmax))
            Nrec = 0.
            Aitrue = []
            sumNrec = 0.
            for irec in range(nbins-1):
                if irec == 0:
                    #Nrec = tree.GetEntries('%s>%.2f && %s<%.2f &&  %s<%.2f' %(self.mtrue, binning[itrue], self.mtrue, binning[itrue+1], self.mrec, binning[irec+1]))
                    Nrec = tree.GetEntries('%s>%.2f && %s<%.2f && %s>%.2f && %s<%.2f' %(self.mtrue, binning[itrue], self.mtrue, binning[itrue+1], self.mrec, binning[irec], self.mrec, binning[irec+1]))
                elif irec == nbins-2:
                    #Nrec = tree.GetEntries('%s>%.2f && %s<%.2f &&  %s>%.2f' %(self.mtrue, binning[itrue], self.mtrue, binning[itrue+1], self.mrec, binning[irec]))
                    Nrec = tree.GetEntries('%s>%.2f && %s<%.2f && %s>%.2f && %s<%.2f' %(self.mtrue, binning[itrue], self.mtrue, binning[itrue+1], self.mrec, binning[irec], self.mrec, binning[irec+1]))
                else:
                    Nrec = tree.GetEntries('%s>%.2f && %s<%.2f && %s>%.2f && %s<%.2f' %(self.mtrue, binning[itrue], self.mtrue, binning[itrue+1], self.mrec, binning[irec], self.mrec, binning[irec+1]))
                Aitrue.append(float(Nrec)/float(Ntrue))
                sumNrec += Nrec
            if self.debug and 0:
                print('itrue, sumNrec/Ntrue, Aitrue =', itrue, sumNrec/Ntrue, Aitrue)
            A.append(Aitrue)
        AmatrixT = np.matrix(A)
        Amatrix = AmatrixT.getT()
        return Amatrix
    def get_isplit_rhomax(self, Vx, htrue):
        n = len(Vx)
        if self.debug:
            print('rhomax n =',n)
            print('Vx =', Vx)
        isplit = -1
        rhomax = -9999.
        for i in range(n):
            dx = htrue.GetBinWidth(i+1)
            if dx < self.res:
                continue
            vii = Vx.item((i,i))
            rhoi = 0.
            for j in range(n):
                vjj = Vx.item((j,j))
                vij = Vx.item((i,j))
                rhoij = vij/math.sqrt(vii*vjj)
                rhoi += rhoij/n
            if rhoi > rhomax and rhoi>0.6:
                rhomax = rhoi
                isplit = i+1
        return isplit

    def get_isplit_ymax(self, hx):
        nbins = hx.GetNbinsX()
        Isplit = -1
        ymax = -9999
        for i in range(nbins):
            y = hx.GetBinContent(i+1)
            dx = hx.GetBinWidth(i+1)
            if y>ymax and dx>=self.res/3.:
                ymax = y
                Isplit = i+1
        return Isplit, ymax
    def get_isplit_dydx(self, hx):
        nbins = hx.GetNbinsX()
        Isplit = -1
        dydxmax = -9999
        for i in range(nbins-1):
            y0 = hx.GetBinContent(i+1)
            y1 = hx.GetBinContent(i+2)
            dx = 0.5*(hx.GetBinWidth(i+1) + hx.GetBinWidth(i+2))
            dydx = (y1-y0)/dx
            if abs(dydx)>dydxmax:
                dydxmax = dydx
                if dydx > 0:
                    Isplit = i+1
                else:
                    Isplit = i+2
        return Isplit
    def get_isplit_peak(self, hx):
        nbins = hx.GetNbinsX()
        Isplit = -1
        ymax = -9999
        for i in range(2, nbins):
            y0 = hx.GetBinContent(i-1)
            y1 = hx.GetBinContent(i)
            y2 = hx.GetBinContent(i+1)
            dx1 = hx.GetBinWidth(i)
            dx0 = 0.5*(hx.GetBinWidth(i-1) + hx.GetBinWidth(i))
            dx2 = 0.5*(hx.GetBinWidth(i) + hx.GetBinWidth(i+1))
            dydx0 = (y1-y0)/dx0
            dydx2 = (y2-y1)/dx2
            if dydx0>0 and dydx2<0 and dx1>self.res/3. and y1> ymax:
                ymax = y1
                Isplit = i
        return Isplit,ymax
    def get_x_fromAb(self, A, b):
        AI = A.getI()
        AIT = AI.getT()
        LVb = []
        n = len(b)
        for i in range(n):
            LVb.append(b.item((i)))
        Vb = np.diag(LVb)
        x = AI.dot(b)
        Vx = AI.dot(Vb.dot(AIT))
        return x, Vx
    def get_b_fromAx(self, A, x):
        return A.dot(x)
    def get_newVx(self, Vx, Isplit, Ddx):
        LnewVx = []
        nbins = len(Vx)
        if self.debug:
            print('nbins =',nbins)
            print('old Vx =', Vx)
            print('Isplit =', Isplit)
            print('Ddx =', Ddx)
        for i in range(nbins):
            v = []
            vii = Vx.item((i,i))
            for j in range(nbins):
                vij = Vx.item((i,j))
                vjj = Vx.item((j,j))
                if i==Isplit-1:
                    if j == Isplit-1:
                        v.append(0.25*vii+Ddx)
                        v.append(0.25*vii-Ddx)
                    else:
                        v.append(0.5*vij)
                else:
                    if j==Isplit-1:
                        v.append(0.5*vij)
                        v.append(0.5*vij)
                    else:
                        v.append(vij)
            LnewVx.append(v)
            if i == Isplit-1:
                v = []
                for j in range(nbins):
                    vij = Vx.item((i,j))
                    if j == Isplit-1:
                        v.append(0.25*vii-Ddx)
                        v.append(0.25*vii+Ddx)
                    else:
                        v.append(0.5*vij)
                LnewVx.append(v)
        newVx = np.matrix(LnewVx)
        if self.debug:
            print('newVx =', newVx)
        return newVx
    def get_newx1(self, x, Isplit, dx):
        if self.debug:
            print('we are in get_b_fromAx')
            print('x =', x, len(x))
        newx = []
        for i in range(len(x)):
            #a = float(x[i][0])
            a = x.item((i))
            if 0:
                print('a =',a)
            if i==Isplit-1:
                newx.append([a*0.5+dx])
                newx.append([a*0.5-dx])
            else:
                newx.append([a])
        newx = np.array(newx)
        if self.debug:
            print('newx =', newx, len(newx))
        return newx
    def get_newx2(self, x, Isplit, dx1, dx2):
        if self.debug and 0:
            print('we are in get_b_fromAx')
            print('x =', x, len(x))
        newx = []
        for i in range(len(x)):
            a = float(x[i][0])
            if 0:
                print('a =',a)
            if i==Isplit-1:
                newx.append([a*0.5+dx1])
                newx.append([a*0.5+dx2])
            else:
                newx.append([a])
        newx = np.array(newx)
        if self.debug and 0:
            print('newx =', newx, len(newx))
        return newx
    def convert2hist(self, hb, x, Vx):
        if self.debug:
            print('hb =', hb)
            print('x =', x)
            print('Vx =', Vx)
        htrue = hb.Clone()
        for i in range(htrue.GetNbinsX()):
            htrue.SetBinContent(i+1, x.item((i)))
            if self.dotoymc:
               htrue.SetBinError(i+1, math.sqrt(Vx.item((i,i))))
            else:
               htrue.SetBinError(i+1, math.sqrt(x.item((i))))
        return htrue
    def get_dlogL(self, A, b, x, hx, isplit, dx=0, eps=0.01):
        logL = self.get_logL(A, b, x, hx, isplit, dx)
        logLp = self.get_logL(A, b, x, hx, isplit, dx+eps)
        dlogL = (logLp - logL)/eps
        return dlogL
    def get_d2logL(self, A, b, x, hx, isplit, dx=0, eps=0.01):
        dlogL = self.get_dlogL(A, b, x, hx, isplit, dx, eps)
        dlogLp = self.get_dlogL(A, b, x, hx, isplit, dx+eps, eps)
        d2logL = (dlogLp - dlogL)/eps
        return d2logL
    def get_logL(self, A, b, x, hx, isplit, dx=0):
        newx = self.get_newx1(x, isplit, dx)
        newb = self.get_b_fromAx(A, newx)
        if self.debug:
            print('we are in get_logL()')
            print('b =',b)
            print('newb =', newb)
        logL = 0.
        nbins = len(b)
        for i in range(nbins):
            if 0:
                print(i, b[i], newb[i])
                print(i, float(b[i]), float(newb[i]))
            db = newb.item((i))-b.item((i))
            if 1.+db/b.item((i)) <= 0:
                print('log(<0)', db, b.item((i)))
                print('newx =', newx)
                print('b =', b)
                print('newb =', newb)
                logL = 0.
                break
            logL += b.item((i))*math.log(1.+db/b.item((i)))-db
        lx = []
        ly = []
        lbins = []
        if isplit==1:
            lbins = [isplit, isplit+1, isplit+2]
        elif isplit==hx.GetNbinsX():
            lbins = [isplit-2, isplit-1, isplit]
        else:
            lbins = [isplit-1, isplit, isplit+1]
        for i in lbins:
            lx.append(hx.GetBinCenter(i))
            ly.append(hx.GetBinContent(i)/hx.GetBinWidth(i))

        local_structure_width = self.get_local_structure_width(lx,ly)
        tau = 1.
        sigma = 0.5*hx.GetBinContent(isplit)
        x1 = 0.
        x4 = 0.
        mu = 0.
        # 1 order regularization
        #if isplit>1:
        #   x1 = hx.GetBinContent(isplit-1)/hx.GetBinWidth(isplit-1)*0.5*hx.GetBinWidth(isplit)
        #if isplit<hx.GetNbinsX():
        #   x4 = hx.GetBinContent(isplit+1)/hx.GetBinWidth(isplit+1)*0.5*hx.GetBinWidth(isplit)
        #mu = (x1-x4)/6.
        chi2 = tau*pow((dx-mu)/sigma,2)
        if local_structure_width > self.res/2.:
           logL -= chi2
        #logL -= chi2
        return logL

    def get_local_structure_width(self, lx=[], ly=[]):
        n = len(lx)
        # solve y = ax^2+bx+c  for 3 points
        # y = y0[1-0.5*(x-mu)^2/sigma^2]
        # sigma^2 = (b^2-4ac)/(8a^2)
        x1 = lx[0]
        x2 = lx[1]
        x3 = lx[2]
        y1 = ly[0]
        y2 = ly[1]
        y3 = ly[2]
        a = - (y1*(x2-x3)+y2*(x3-x1)+y3*(x1-x2))/((x1-x2)*(x2-x3)*(x3-x1))
        b = (y1-y2)/(x1-x2) - a*(x1+x2)
        c = y1- a*x1*x1-b*x1
        sigmasq = 9999
        if a!= 0:
           sigmasq = (b*b-4*a*c)/(8*a*a)
        #mu = b/(-2.*a)
        sigma = math.sqrt(abs(sigmasq))
        sign = sigmasq/abs(sigmasq)
        #return sigma*sign
        return sigma
    def unfold_core(self, sigtree, datatree, nsplit=-1, dotoymc=0):
        if self.debug:
            print('signal tree', sigtree, sigtree.GetEntries())
            print('data tree', datatree, datatree.GetEntries())
        if nsplit>=0:
            self.nsplit = nsplit
        binning = self.initbinning
        hb, b = self.get_b(datatree, binning, dotoymc)
        A = self.get_response_matrix(sigtree, binning)
        x, Vx = self.get_x_fromAb(A, b)
        x0 = x
        if not dotoymc:
            self.Lisplit = []
        for i in range(self.nsplit):
            if dotoymc and i>=len(self.Lisplit):
                print('this toymc over!')
                break
            if self.debug or 0:
                print('start split',i,'x =',x)
            hx = self.convert2hist(hb, x, Vx)
            if self.debug or 0:
                sumVxii = 0.
                for a in range(len(Vx)):
                    for b in range(len(Vx)):
                        sumVxii += Vx.item((a,b))
                print('sumVxii =', sumVxii)
            isplit = -1
            if dotoymc:
                isplit = self.Lisplit[i]
            else:
                isplit,ymax = self.get_isplit_peak(hx)
                if isplit < 0:
                    isplit,ymax = self.get_isplit_ymax(hx)
                if isplit<0:
                    print('HAHA, no isplit found!!!!!!!!!!!!!')
                    break
                self.Lisplit.append(isplit)
            binning = self.get_binning(binning, isplit)
            #self.binning = binning
            hb, b = self.get_b(datatree, binning, dotoymc)
            A = self.get_response_matrix(sigtree, binning)
            Nnewtons = 20
            dx = 0 # initial value for dx in the Newton method
            bestfit_d2logL = 0.
            for itr in range(Nnewtons):
                dx0 = dx
                dlogL = self.get_dlogL(A, b, x, hx, isplit, dx0, 0.01)
                d2logL = self.get_d2logL(A, b, x, hx, isplit, dx0, 0.01)
                if d2logL < 0:
                    dx = dx0 - dlogL/d2logL
                    bestfit_d2logL = d2logL
                elif d2logL > 0:
                    if self.debug:
                        print('d2logL > 0, it is not maximizing logL at',itr+1,'iteraction')
                    break
                if abs(dx-dx0)<0.001*dx0:
                    if self.debug:
                        print('finish Newton method with |dx-dx0|<0.001*dx0 at',itr+1,'iteration')
                    break
                if abs(dx-dx0)<0.01:
                    if self.debug:
                        print('finish Newton method with |dx-dx0|<0.01 at', itr+1,'iteration')
                    break
            if abs(dx) > 0.5*x.item((isplit-1)):
                #print('warning: dx>0.5*x, set dx to its limit')
                #dx = 0.5*x.item((isplit-1))*dx/abs(dx)
                print('warning: dx>0.5*x, set dx to 0')
                dx = 0.
            # now dx is the best dx, we need to update x and Vx
            x = self.get_newx1(x, isplit, dx)
            Ddx = 1./(-bestfit_d2logL)
            #Ddx = 0.25*Vx.item((isplit-1, isplit-1))
            Vx= self.get_newVx(Vx, isplit, Ddx)
            if self.debug:
                print('end split',i,'x =',x)
        self.Vx = Vx
        # get the meas. distribution with the final binning
        hb, b = self.get_b(datatree, binning, dotoymc)
        # prepare histograms for final output
        #hrec = hb.Clone('hrec_clone')
        #htrue = self.get_htrue(sigtree, binning)
        hout = self.convert2hist(hb, x, Vx)
        #binning = binning
        return  x, hout, binning

    def do_unfold(self, nsplit = -1):
        sigfile = TFile(self.sig_filepath, 'read')
        sigtree = sigfile.Get(self.sig_treename)
        datafile = TFile(self.data_filepath, 'read')
        datatree = datafile.Get(self.data_treename)
        if self.debug:
            print('tree, nentries =', tree.GetName(), tree.GetEntries())
        x, self.hout, self.binning = self.unfold_core(sigtree, datatree, nsplit)
        self.hrec, b = self.get_b(datatree, self.binning, 0)
        self.htrue = self.get_htrue(sigtree, self.binning)

        if self.dotoymc:
            self.do_toymc(sigtree, self.hrec, self.hout)
            self.hout = self.convert2hist(self.hrec, x, self.Vx)

        self.Lrec = self.convert2list(self.hrec)
        self.Ltrue = self.convert2list(self.htrue)
        self.Lout = self.convert2list(self.hout)
        sigfile.Close()
        datafile.Close()
        # calculate covariance matrix and correlation matrix
        self.cal_Rx(self.Vx)
    def do_toymc(self, sigtree, hrec, hout0):
        print('We are doing toy mc...')
        if self.ntoys < 0:
            print('WARNING: ntoys < 0, set it to 1')
            self.ntoys = 1
        nevents = int(hrec.Integral())
        list_hout = []
        toymcfile = TFile('toymcfile.root', 'recreate')
        for itoy in range(self.ntoys):
            print('doing toymc ', itoy, '...')
            datatree = TTree(self.data_treename, 'toymc_'+str(itoy))
            mrec = array('f', [0])
            datatree.Branch('mrec', mrec, 'mrec/F')
            for i in range(nevents):
                mrec[0] = hrec.GetRandom()
                datatree.Fill()
            datatree.Write()
            x, hout, binning = self.unfold_core(sigtree, datatree, self.nsplit, dotoymc=1)
            list_hout.append(hout)
        if self.debug or 0:
            self.showhist(hout0)
            for hout in list_hout:
                print('toymc results', hout)
                self.showhist(hout)
        self.Vx = self.get_newVx_from_toymc(hout0, list_hout)
        toymcfile.Close()
    def get_newVx_from_toymc(self, hout0, list_hout):
        nbins = hout0.GetNbinsX()
        ntoys = len(list_hout)
        list_avg_x = []
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>toy mc details: checking bias ...')
        for i in range(nbins):
            avg_x = 0.
            for hout in list_hout:
                avg_x += hout.GetBinContent(i+1)
            avg_x = float(avg_x)/ntoys
            print(i, hout0.GetBinContent(i+1), avg_x)
            list_avg_x.append(avg_x)
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>toy mc details: checking correlation ...')
        Vx = []
        for i in range(nbins):
            covi = []
            for j in range(nbins):
                covij = 0.
                for hout in list_hout:
                    xi = hout.GetBinContent(i+1)
                    xj = hout.GetBinContent(j+1)
                    covij += xi*xj
                covij = float(covij)/ntoys - list_avg_x[i]*list_avg_x[j]
                covi.append(covij)
                print(i,j, covij)
            Vx.append(covi)
        Vx = np.matrix(Vx)
        return Vx

    def get_output(self):
        return self.binning, self.Lrec, self.Ltrue, self.Lout
    def get_output_hist(self):
        return self.binning, self.hrec, self.htrue, self.hout
    def cal_Rx(self, Vx):
        n = len(Vx)
        hVx = TH2F('hVx', '', n, 0, n, n, 0, n)
        hRx = TH2F('hRx', '', n, 0, n, n, 0, n)
        #print('1, Vx =', Vx)
        Rx = np.zeros((n,n))
        #print('2, Vx =', Vx)
        for i in range(n):
            for j in range(n):
                vij = Vx.item((i,j))
                vii = Vx.item((i,i))
                vjj = Vx.item((j,j))
                rij = vij/math.sqrt(vii*vjj)
                #Rx.itemset((i,j), rij)
                Rx[i,j] = rij
                hVx.SetBinContent(i+1, j+1, vij)
                hRx.SetBinContent(i+1, j+1, rij)
        #print('after Vx =', Vx)
        #print('Rx =', Rx)
        self.Rx = Rx
        self.hVx = hVx
        self.hRx = hRx
    def get_Rx(self):
        #return self.hVx, self.Vx
        return self.hRx, self.Rx
    def convert2list(self,h):
        L = []
        for i in range(h.GetNbinsX()):
            L.append([h.GetBinContent(i+1), h.GetBinError(i+1)])
        return L
    def showhist(self, h):
        print('showhist', h)
        nbins = h.GetNbinsX()
        print(h, nbins)
        for i in range(nbins):
            print(i, h.GetBinContent(i+1), '+/-', h.GetBinError(i+1))


