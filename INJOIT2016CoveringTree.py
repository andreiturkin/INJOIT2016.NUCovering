import abc
import numpy as np
from math import sqrt
from ete3 import Tree
#Plotting
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#Images
from PIL import Image
# Date and Time
import datetime

# =============================================================================
# Extension modules
# =============================================================================
import pyOpt

class Rect:
    def __init__(self, left, top, width, height):
        self.left = left
        self.right = left + width
        self.top = top
        self.bottom = top + height
        self.width = width
        self.height = height
        self.centerx = left + width/2
        self.centery = top + height/2
        self.center = (left + width/2, top + height/2)
    
    def __str__(self):
        return '<Rect: {}, {}, {}, {}>'.format(self.left, self.top, self.width, self.height)  

class CoveringTree:
############################################################################################
# Constructor
############################################################################################    
    __metaclass__ = abc.ABCMeta
    def __init__(self, l0, l1_bounds, l2_bounds, idelta=0):
        #Initialize Base
        self.l0 = l0
        #Initialize the Bounds
        self.__l1_bounds = l1_bounds
        self.__l2_bounds = l2_bounds
        
        #Define Initial Rectangle P
        left = -self.__l1_bounds[1]
        top = 0
        width = self.__l1_bounds[1]+self.l0+self.__l2_bounds[1]
        height = min(self.__l1_bounds[1],self.__l2_bounds[1])
        
        #Initialize initial Space where the workspace lie
        self.__Xspace = Rect(left, top, width, height)
        #Initialize the Root
        self.__initTree(self.__Xspace)
        #Initialize the minimal size of the rectangle
        self.__delta = idelta
        
        #Initialize plotting facilities
        self.__fig = plt.figure()
        
        self.__ax = self.__fig.add_subplot(111)
        self.__ax.axis('scaled')
        self.__ax.axis([self.__Xspace.left, self.__Xspace.right, self.__Xspace.top, self.__Xspace.bottom ])
        
        self.__tleveltext = self.__ax.text(0.98,0.9, 'Tree Level = {}'.format(0),\
                       verticalalignment='center', \
                       horizontalalignment='right', \
                       fontsize=13,\
                       transform = self.__ax.transAxes)
        
        self.__curdiam = self.__ax.text(0.02,0.9, 'd(Rectangle) = {}'.format(round(self.__d(self.__Xspace),4)),\
                       verticalalignment='center', \
                       horizontalalignment='left', \
                       fontsize=13,\
                       transform = self.__ax.transAxes)

    @abc.abstractmethod
    def getMaxVal(self, xbounds, ybounds, diam):
        raise NotImplementedError
    @abc.abstractmethod
    def getMinVal(self, xbounds, ybounds, diam):
        raise NotImplementedError
############################################################################################
# Private Members
############################################################################################        
    def __vSplitter(self, iRect):
        newleft1 = iRect.left
        newtop1 = iRect.top
        newwidth1 = iRect.width/2.0
        newheight1 = iRect.height
        Rleft = Rect(newleft1, newtop1, newwidth1, newheight1)
        
        newleft2 = iRect.left + iRect.width/2.0
        newtop2 = iRect.top
        newwidth2 = iRect.width/2.0
        newheight2 = iRect.height
        Rright = Rect(newleft2, newtop2, newwidth2, newheight2)
        return Rleft, Rright
    
    def __hSplitter(self, iRect):
        newleft1 = iRect.left
        newtop1 = iRect.top
        newwidth1 = iRect.width
        newheight1 = iRect.height/2.0
        Rleft = Rect(newleft1, newtop1, newwidth1, newheight1)
        
        newleft2 = iRect.left
        newtop2 = iRect.top + iRect.height/2.0
        newwidth2 = iRect.width
        newheight2 = iRect.height/2.0
        Rright = Rect(newleft2, newtop2, newwidth2, newheight2)
        return Rleft, Rright
    
    def __d(self, iRect):
        return sqrt(iRect.width**2.0 + iRect.height**2.0)
    
    def g1(self, x):
        return x[0]**2.0 + x[1]**2.0 - (self.__l1_bounds[1]**2.0) 
    
    def g2(self, x):
        return self.__l1_bounds[0]**2.0 - (x[0]**2.0) - (x[1]**2.0)
    
    def g3(self, x):
        return ((x[0] - self.l0)**2.0) + (x[1]**2.0) - (self.__l2_bounds[1]**2.0)
         
    def g4(self, x):
        return self.__l2_bounds[0]**2.0 - ((x[0] - self.l0)**2.0) - (x[1]**2.0)
    
    def g3m(self, x):
        return np.array([(x[0]**2.0) + (x[1]**2.0) - (self.__l2_bounds[1]**2.0)])
        
    def g4m(self, x):
        return np.array([self.__l2_bounds[0]**2.0 - (x[0]**2.0) - (x[1]**2.0)])
    
    def phi(self, x):
        return max(self.g1(x), self.g2(x), self.g3(x), self.g4(x))
    
    def __analyseRect(self, iRect):
        
        xmin = iRect.left
        xmax = iRect.left + iRect.width
        ymin = iRect.top
        ymax = iRect.top + iRect.height
        
        maxval = self.getMaxVal((xmin,xmax), (ymin,ymax), self.__d(iRect)) 
        #The whole rectangle is a part of the solution -> save it
        if (maxval < 0):
            #mark it as in range
            inrange = True
            return False, inrange 
                
        minval = self.getMinVal((xmin,xmax), (ymin,ymax), self.__d(iRect)) 
        #There is no solution for the rectangle -> get rid of it
        if (minval > 0):
            #mark it as out of range
            inrange = False
            return False, inrange
        
        #The rectangle should be processed further
        return True, False

    def __addToTree(self, motherNode, iRect1, iRect2, childNodeLevel):
        # and add the nodes as children.
        oNode2 = motherNode.add_child(name='{}'.format(childNodeLevel))
        oNode1 = motherNode.add_child(name='{}'.format(childNodeLevel))
        #add features
        oNode2.add_feature('Rect',iRect2)
        oNode1.add_feature('Rect',iRect1)
        
    def __getNewRect(self, iRect, level):        
        (oRleft,oRright) = self.__vSplitter(iRect) if (level%2==0) else self.__hSplitter(iRect)
        return (oRleft,oRright)
    
    def __initTree(self, Xspace):
        self.__sTree = Tree('0;') #name here is the level of the tree
        motherNode = self.__sTree.search_nodes(name='0')[0]
        motherNode.add_feature('Rect',Xspace)
    
    def __drawRect(self, iRect, fillIt, PlotEdges=True, inQI=False, inQE=True):
        if(PlotEdges):
            #Internal
            if inQI and inQE:
                edgeColor = 'black'
                LineStyle='solid'
                LineWidth = 1
                Alpha=0.3
            #External
            if inQE and (not inQI):
                edgeColor = 'red'
                LineStyle='solid'
                LineWidth = 1
                Alpha=None
            #Out of range
            if (not inQE) and (not inQI):
                edgeColor = 'green'
                LineStyle='solid'
                LineWidth = 1
                Alpha=None
            
            self.__ax.add_patch(
                          patches.Rectangle(
                                            (iRect.left, iRect.top),   # (x,y)
                                            iRect.width,          # width    
                                            iRect.height,         # height
                                            fill = inQI,
                                            alpha = Alpha,
                                            linestyle = LineStyle,
                                            edgecolor = edgeColor,  
                                            lw = LineWidth)
                          )
        else:
            self.__ax.add_patch(
                          patches.Rectangle(
                                            (iRect.left, iRect.top),   # (x,y)
                                            iRect.width,          # width    
                                            iRect.height,         # height
                                            fill = fillIt,
                                            edgecolor = 'none')
                          )
        plt.draw()
        
############################################################################################
# Public Members
############################################################################################    
    def getCovering(self, maxLevels, saveasmovie=True):
        
        cdRect = self.__d(self.__Xspace)
        print 'The diameter of the initial rectangle is {}\n'.format(cdRect)

        bExit = False
        for curLevel in range(0, maxLevels):
            print 'Processing level {}'.format(curLevel)
            #pause(0.000001)
            #Get all the rectangles that are on some level of the tree
            curLevelNodes = self.__sTree.get_leaves_by_name(name='{}'.format(curLevel))
            #Loop over the rectangles
            for curLevelNode in curLevelNodes:
                #Get a rectangle from the tree level
                oRect = curLevelNode.Rect
                #Save current rectangle diameter
                if self.__d(oRect) < cdRect:
                    cdRect = self.__d(oRect)
                    print 'Current level diameter of the rectangle is {}\n'.format(cdRect)
                
                inQE = False
                inQI = False
                #The diameter of the rectangle is less than or equal to the predefined delta value       
                #see eq. 2.6: d(P^(i)) <= \delta
                if self.__d(oRect) <= self.__delta:
                    #It is too small to decide upon -> save it as if it was in range
                    cont = False
                    inrange = True
                    inQE = True
                    inQI = False
                    #Return the result on the next iteration 
                    bExit = True
                #Otherwise
                else:
                    #Analyze it
                    #see eq. 2.4 and 2.5
                    (cont, inrange) = self.__analyseRect(oRect)
                    if inrange: 
                        inQI = True
                        inQE = True
                #Save the obtained results
                if cont and (curLevel < maxLevels-1):  
                    (oRleft,oRright) = self.__getNewRect(oRect,curLevel)
                    self.__addToTree(curLevelNode, oRleft, oRright, curLevel + 1)
                else:
                    #save results to the analyzed node
                    curLevelNode.add_feature('Inrange',inrange)
                    curLevelNode.add_feature('inQI',inQI)
                    curLevelNode.add_feature('inQE',inQE)
                    
            #All of the rectangles could be obtained on the next iterations are too small
            #so break it
            if bExit:
                print 'The result is obtained for {} levels'.format(curLevel)
                break
        
        #plt.show()

    def saveCoveringAsImage(self, fileName='./Images/{0}__{1:02d}_{2:02d}_{3:02d}_covering.jpeg'.format(datetime.date.today(), \
                                                           datetime.datetime.now().hour,\
                                                           datetime.datetime.now().minute,\
                                                           datetime.datetime.now().second),\
                                                           ResOnly = False, Grayscale = False):
        plt.cla()
        
        for leaf in self.__sTree.iter_leaves():
            if(ResOnly):
                #Draw the rectangle without edges
                self.__drawRect(leaf.Rect, leaf.Inrange, False)
            else:
                #Draw the rectangle with edges
                self.__drawRect(leaf.Rect, leaf.Inrange, True, leaf.inQI, leaf.inQE)
            
        plt.draw()
        plt.pause(1)
        
        if (Grayscale):
            self.__fig.savefig('./Images/temp.png', dpi = 600)
            Image.open('./Images/temp.png').convert("L").save(fileName)
        else:
            self.__fig.savefig(fileName, dpi = 600)


class CoveringTreeGlobOpt(CoveringTree):
    def __init__(self, l0, l1_bounds, l2_bounds, idelta=0):
        CoveringTree.__init__(self, l0, l1_bounds, l2_bounds, idelta)
        
    def objfunc(self, x):
        f = max(self.g1(x), self.g2(x), self.g3(x), self.g4(x))
        g = []

        fail = 0
        return f,g, fail

    def getMinVal(self, xbounds, ybounds, diam):
        xmin = xbounds[0]
        xmax = xbounds[1]
        ymin = ybounds[0]
        ymax = ybounds[1]
        
        line = 'min(phi(x)), {}<=x1<={}, {}<=x2<={}'.format(xmin,xmax,ymin,ymax)
        opt_prob = pyOpt.Optimization(line, self.objfunc)
        opt_prob.addObj('phi')
        opt_prob.addVar('x1','c',lower=xmin,upper=xmax,value=(xmax+xmin)/2.0)
        opt_prob.addVar('x2','c',lower=ymin,upper=ymax,value=(ymax+ymin)/2.0)
         
        nsga2 = pyOpt.MIDACO()
        nsga2.setOption('IPRINT',-1)
        nsga2(opt_prob)
        
        #MINPHI(x) = MINMAX(g1(x),g2(x),g3(x),g4(x)
        return opt_prob.solution(0)._objectives[0].value

    def getMaxVal(self, xbounds, ybounds, diam):
        xmin = xbounds[0]
        xmax = xbounds[1]
        ymin = ybounds[0]
        ymax = ybounds[1]
         
        #MAX
        #g1(x1,x2)
        g1a1max = max(abs(xmin),abs(xmax))
        g1a2max = max(abs(ymin),abs(ymax))
        #g2(x1,x2)
        g2a1max = min(abs(xmin),abs(xmax))
        g2a2max = min(abs(ymin),abs(ymax))
        #g3(x1,x2)
        g3a1max = max(abs(xmin-self.l0),abs(xmax-self.l0))
        g3a2max = max(abs(ymin),abs(ymax))
        #g4(x1,x2)
        g4a1max = min(abs(xmin-self.l0),abs(xmax-self.l0))
        g4a2max = min(abs(ymin),abs(ymax))
        
        #MAXPHI(x) = MAXMAX(g1(x),g2(x),g3(x),g4(x) 
        return max(self.g1((g1a1max,g1a2max)),self.g2((g2a1max,g2a2max)),\
                   self.g3m((g3a1max,g3a2max)),self.g4m((g4a1max,g4a2max)))     

class CoveringTreeAppx(CoveringTree):

    def getMinVal(self, xbounds, ybounds, diam):
        xmin = xbounds[0]
        xmax = xbounds[1]
        ymin = ybounds[0]
        ymax = ybounds[1]
        
        #L = sup_{x in P}(||grad(phi(x))||)=sup_{x in P}(||grad(gi(x))||)
        half_L = max(sqrt(max(abs(xmin),abs(xmax))**2 + max(abs(ymin),abs(ymax))**2),\
                     sqrt(max(abs(xmin-self.l0), abs(xmax-self.l0))**2 + max(abs(ymin),abs(ymax))**2)) 
        
        #MINPHI(x) = MINMAX(g1(x),g2(x),g3(x),g4(x)
        return self.phi(((xmin+xmax)/2,(ymin+ymax)/2))-half_L*diam
    
    def getMaxVal(self, xbounds, ybounds, diam):
        xmin = xbounds[0]
        xmax = xbounds[1]
        ymin = ybounds[0]
        ymax = ybounds[1]
        
        #L = sup_{x in P}(||grad(phi(x))||)=sup_{x in P}(||grad(gi(x))||)
        half_L = max(sqrt(max(abs(xmin),abs(xmax))**2 + max(abs(ymin),abs(ymax))**2),\
                     sqrt(max(abs(xmin-self.l0), abs(xmax-self.l0))**2 + max(abs(ymin),abs(ymax))**2)) 
        
        #MINPHI(x) = MINMAX(g1(x),g2(x),g3(x),g4(x)
        return self.phi(((xmin+xmax)/2,(ymin+ymax)/2)) + half_L*diam

