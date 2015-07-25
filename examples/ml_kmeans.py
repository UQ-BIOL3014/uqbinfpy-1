def test_kmeans():
    import matplotlib.pyplot as plt # plotting/visualization
    from mpl_toolkits.mplot3d import Axes3D

    data  = numpy.empty((120,3))
    data[ 0: 40,:] = numpy.random.randn(40,3) * 0.53 + 1 # red
    data[40: 80,:] = numpy.random.randn(40,3) * 0.20 - 1 # green
    data[80:120,:] = numpy.random.randn(40,3) * 0.57 + 2 # blue

    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax1 = fig.add_subplot(1,2,1,projection='3d')
    ax1.plot(data[ 0: 40,0],data[ 0: 40,1],data[ 0: 40,2],'r.') # red
    ax1.plot(data[40: 80,0],data[40: 80,1],data[40: 80,2],'g.') # green
    ax1.plot(data[80:120,0],data[80:120,1],data[80:120,2],'b.') # blue

    km = KMeans(data)
    km.train(3)
    ax2 = fig.add_subplot(1,2,2,projection='3d')
    print "Means:"
    for i in range(3):
        print km.means[i]
    ax2.plot([km.means[0,0]],[km.means[0,1]],[km.means[0,2]],'ro') # red
    ax2.plot([km.means[1,0]],[km.means[1,1]],[km.means[1,2]],'go') # green
    ax2.plot([km.means[2,0]],[km.means[2,1]],[km.means[2,2]],'bo') # blue

    print "Classifications:"
    for i in range(len(data)):
        print data[i], km.classify(data[i])
        if km.classify(data[i]) == 0:
            ax2.plot([data[i,0]],[data[i,1]],[data[i,2]],'r.') # red
        elif km.classify(data[i]) == 1:
            ax2.plot([data[i,0]],[data[i,1]],[data[i,2]],'g.') # green
        else:
            ax2.plot([data[i,0]],[data[i,1]],[data[i,2]],'b.') # blue
    print "Close plot to exit."
    plt.show()

if __name__=='__main__':
    test_kmeans()



