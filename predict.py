def predict(gene, model, visualise=True, verbose=1):
    """
    Predict new gene (this is NOT a method)
    """
    probs = model.predict(gene.data, batch_size=None)
    if len(gene.targets.shape) == 1:
        probs = [:,0]
        MAP = np.where(probs < 0.5, 0., 1.)
        MLE = np.average(gene.classes, weights = probs[0])

        MAP = self.gene.classes[np.argmax(probs)]
        MLE = np.average(self.gene.classes, weights = probs[0])
        BF = (1 - probs[0][0]) / probs[0][0]
        samples_distr = np.random.choice(self.gene.classes, size = 100000, replace = True, p = probs[0])
        HPD = pymc3.stats.hpd(samples_distr, alpha = 0.05)
        if verbose > 0:
            print('MLE (', str(round(MLE,2)), ')')
            print('MAP (', str(MAP), ')')
            print('BF (', str(round(BF,2)), ')')
            print('HPD (', str(HPD), ')')

        if visualise == True:
            plt.figure()
            tick_marks = self.gene.classes
            cen_tick = self.gene.classes
            plt.hist(samples_distr, color='#a6bddb', bins=len(self.gene.classes), density=True)
            plt.xlim([self.gene.classes.min(), self.gene.classes.max()])
            plt.xticks(cen_tick, cen_tick, rotation=45, fontsize=10)
            plt.yticks(fontsize=10)
            plt.ylabel('Probability', fontsize=12)
            plt.xlabel('Parameter', fontsize=12)
            plt.title('Sampled posterior distribution')
            plt.grid(True)
            plt.show()
        return (probs, (MAP, MLE, BF, HPD))


