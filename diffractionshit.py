   def fourier(M, r):
        Mx = M[0]
        My = M[1]
        Mz = M[2]
        print np.amin(Mz), np.amax(Mz)
        colors = [plt.cm.viridis(x) for x in 0.5+Mz/(2*np.amax(Mz))]
        plt.subplot()
        ax1 = plt.subplot(121)
        # Z = np.transpose(np.transpose(colors)[:-1])
        # im = Z.reshape((len(xvals), len(yvals), 3))
        # ax1.imshow(gaussian_filter(np.rot90(im), sigma = 1), interpolation = 'bilinear',
        # alpha = 0.9, origin = 'upper')
        # ax1.quiver(r[0], -r[1]+100, Mx, My)
        # plt.xticks([])
        # plt.yticks([])

        ax1 = plt.subplot(121)
        plt.xticks([])
        plt.yticks([])
        img = Mz.reshape((len(xvals), len(yvals)))
        f = np.fft.fft2(img)
        f = np.fft.fftshift(f)
        magnitude_spectrum = np.absolute(f)

        ax1.imshow(np.rot90(magnitude_spectrum))  # , interpolation='hamming')

        img = (Mz*window).reshape((len(xvals), len(yvals)))
        ax2 = plt.subplot(122)
        plt.xticks([])
        plt.yticks([])
        f = np.fft.fft2(img)
        f = np.fft.fftshift(f)
        magnitude_spectrum = np.absolute(f)
        ax2.imshow(np.rot90(magnitude_spectrum), interpolation='hamming')
        plt.show()

    def formfactor(m, r):
        M = np.zeros((3, len(xvals), len(yvals)), dtype=np.complex128)
        for i in xrange(len(m)):
            M[i] = np.fft.fftshift(np.fft.fft2(
                (m[i]*window).reshape((len(xvals), len(yvals)))))
        ki = np.array([0, 1, 0], dtype=np.complex128)
        kf = np.array([0, 1, 0], dtype=np.complex128)  # defines transmission
        M = M.transpose()
        Isig = np.zeros((len(yvals), len(xvals)))
        Ipi = np.zeros((len(yvals), len(xvals)))
        IDic = np.zeros((len(yvals), len(xvals)))

        for i in range(len(M)):
            for j in range(len(M[i])):
                # np.sqrt((i-len(xvals)/2.)**2+(j-len(yvals)/2.)**2) # x and y might be other way round...
                q = 1e-7*np.array([i-len(xvals)/2., (j-len(yvals)/2.), 0])
                kf = q + kf
                kf = kf / np.linalg.norm(kf)
                Isig[i][j] = np.absolute(np.vdot(kf, M[i][j]))
                Ipi[i][j] = np.absolute(
                    np.vdot(ki, M[i][j]) + np.vdot(np.cross(ki, kf), M[i][j]))
                IDic[i][j] = np.real(
                    j*np.vdot(np.vdot(kf, np.conj(M[i][j]))*np.cross(ki, kf), M[i][j]))
        plt.figure()
        plt.subplot(131)
        plt.title('$\sigma$')
        plt.xticks([])
        plt.yticks([])
        plt.imshow((Isig))
        plt.subplot(132)
        plt.title('$\pi$')
        plt.xticks([])
        plt.yticks([])
        plt.imshow((Ipi))
        plt.subplot(133)
        plt.title('$Dichroism$')
        plt.xticks([])
        plt.yticks([])
        plt.imshow((IDic))
        plt.show()
