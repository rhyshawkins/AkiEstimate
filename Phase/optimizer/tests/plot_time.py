
import glob
import matplotlib.pyplot as P


if __name__ == '__main__':

    stationpair = 'HOT05_HOT25'
    files = glob.glob('Final_%s_*/opt.time' % stationpair)

    x = []
    y = []
    
    for fn in files:
        f = open(fn, 'r')
        lines = f.readlines()
        f.close()

        seconds = float(lines[0].split()[0][:-4])
        label = int(fn.split('/')[0][len('Final_%s_' % stationpair):])
        

        x.append(label)
        y.append(seconds)

    z = zip(x, y)
    z.sort()
    
    sx, sy = zip(*z)

    fig, ax = P.subplots()
    
    ax.plot(sx, sy)

    P.show()
    


        
