from INJOIT2016CoveringTree import CoveringTreeGlobOpt
from INJOIT2016CoveringTree import CoveringTreeAppx
import timeit
import datetime

def TepeeExp():
    deltas = [0.5, 0.06]

    # Kite like workspace
    for iDelta in deltas:
        print 'Delta = {}'.format(iDelta)
        TepeePlotting(iDelta, False) 

def TepeePlotting(iDelta, ShowResOnly=False):
    
    cTree = CoveringTreeAppx(4, [2, 12], [4, 12], iDelta)

    maxLevels = 64   
    t = timeit.Timer(lambda: cTree.getCovering(maxLevels,False))
    exectime = t.timeit(number=1)
    print 'Execution time: {}'.format(exectime)
    
    cTree.saveCoveringAsImage('./Images/Tepee_{0}__{1:02d}_{2:02d}_{3:02d}_covering_{4}.jpeg'.format(datetime.date.today(), \
                                                           datetime.datetime.now().hour,\
                                                           datetime.datetime.now().minute,\
                                                           datetime.datetime.now().second,\
                                                           iDelta))
    cTree = CoveringTreeGlobOpt(4, [2, 12], [4, 12], iDelta)
  
    t = timeit.Timer(lambda: cTree.getCovering(maxLevels,False))
    exectime = t.timeit(number=1)
    print 'Execution time: {}'.format(exectime)
     
    cTree.saveCoveringAsImage('./Images/Tepee_{0}__{1:02d}_{2:02d}_{3:02d}_covering_{4}.jpeg'.format(datetime.date.today(), \
                                                           datetime.datetime.now().hour,\
                                                           datetime.datetime.now().minute,\
                                                           datetime.datetime.now().second,\
                                                           iDelta))


if __name__ == '__main__':
    TepeeExp()

