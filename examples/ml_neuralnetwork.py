if __name__ == '__main__0':
    im, om = readDenseDataFile('xor.trn')
    nrows, ninputs = im.shape
    nrows, noutputs = om.shape
    nn = NN(ninputs, 3, noutputs)
    for i in range(100):
        print nn.train(im, om, 0.1, 100)

