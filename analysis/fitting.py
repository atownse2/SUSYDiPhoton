import os

import ROOT

import numpy as np
import pandas as pd

top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
plot_dir = f"{top_dir}/outputs/plots/egtf/fits/"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

from analysis.utils.logger import Logger
l = Logger()


class Linear:
    def __init__(self, x):
        self.intercept = ROOT.RooRealVar("c0", "Y intercept of linear bkg", 1, 0, 500)
        self.slope = ROOT.RooRealVar("c1", "Slope of linear bkg", 0.1, 0, 10)
        self.pdf = ROOT.RooGenericPdf("linear", "linear background", "c0 + c1*x", ROOT.RooArgList(self.intercept,self.slope,x))

class Voigt:

    def __init__(self, x, data, min_val=1e-1, max_val=1e7):
        self.name = 'Voigt'

        mode = np.mean(data)
        std = np.std(data)

        self.mean = ROOT.RooRealVar('mean', 'mean', mode, 0.5*mode, 1.5*mode)
        self.width = ROOT.RooRealVar('width', 'width', 0.8*std, 0.2*std, 1.3*std)
        self.sigma = ROOT.RooRealVar('sigma', 'sigma', 0.8*std, 0.2*std, 1.3*std)
        self.pdf = ROOT.RooVoigtian('voigt', 'voigt', x, self.mean, self.width, self.sigma)

    def legend(self, legend):
        legend.AddEntry(0, f'mean: {self.mean.getVal():.2f}', '')
        legend.AddEntry(0, f'width: {self.width.getVal():.2f}', '')
        legend.AddEntry(0, f'sigma: {self.sigma.getVal():.2f}', '')

class CrystalBall:
    def __init__(self, x, data, min_val=1e-1, max_val=1e7):
        self.name = 'CrystalBall'

        mode = np.mean(data)
        std = np.std(data)

        self.mean = ROOT.RooRealVar('mean', 'mean', mode, 0.5*mode, 1.5*mode)
        self.sigma = ROOT.RooRealVar('sigma', 'sigma', 0.8*std, 0.2*std, 1.3*std)
        self.alphaL = ROOT.RooRealVar('alphaL', 'alphaL', 1, min_val, 50)
        self.nL = ROOT.RooRealVar('nL', 'nL', 2, min_val, 100)
        self.alphaR = ROOT.RooRealVar('alphaR', 'alphaR', 1, min_val, 50)
        self.nR = ROOT.RooRealVar('nR', 'nR', 2, min_val, 100)
        self.pdf = ROOT.RooCrystalBall('cb', 'cb', x, self.mean, self.sigma, self.alphaL, self.nL, self.alphaR, self.nR)

    def legend(self, legend):
        legend.AddEntry(0, f'mean: {self.mean.getVal():.2f}', '')
        legend.AddEntry(0, f'sigma_range: {self.sigma.getMin():.2f} - {self.sigma.getMax():.2f}', '')
        legend.AddEntry(0, f'sigma: {self.sigma.getVal():.2f}', '')
        legend.AddEntry(0, f'alphaL: {self.alphaL.getVal():.2f}', '')
        legend.AddEntry(0, f'nL: {self.nL.getVal():.2f}', '')
        legend.AddEntry(0, f'alphaR: {self.alphaR.getVal():.2f}', '')
        legend.AddEntry(0, f'nR: {self.nR.getVal():.2f}', '')

def fit(
    df : pd.DataFrame,
    name : str,
    fit_type : str,
    signal_model = CrystalBall,
    background_model = Linear,
    plot = True,
):
    assert fit_type in ['binned', 'unbinned']
    
    vals = df['invariant_mass'].to_numpy()
    weights = df['weight'].to_numpy()

    x = ROOT.RooRealVar("x", "Invariant Mass", np.floor(vals.min()), np.ceil(vals.max()))

    bkg = background_model(x)
    sig = signal_model(x, vals)

    nsig = ROOT.RooRealVar("nsig", "Number of signal events", 0, sum(weights))
    nbkg = ROOT.RooRealVar("nbkg", "Number of background events", 0, sum(weights))

    model = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(sig.pdf, bkg.pdf), ROOT.RooArgList(nsig,nbkg))

    if fit_type == 'binned':
        histogram = ROOT.TH1D(name, name, 32, 50, 130)
        histogram.Sumw2()
        for val, weight in zip(vals, weights):
            histogram.Fill(val, weight)
        data = ROOT.RooDataHist("data", "data", ROOT.RooArgList(x), histogram)
    elif fit_type == 'unbinned':
        tree = ROOT.TTree('tree', 'tree')
        v = np.zeros(1, dtype=float)
        w = np.zeros(1, dtype=float)
        tree.Branch('v', v, 'v/D')
        tree.Branch('weights', w, 'weights/D')
        for val, weight in zip(vals, weights):
            v[0] = val
            w[0] = weight
            tree.Fill()

        data = ROOT.RooDataSet('data', 'data', tree, ROOT.RooArgSet(x), 'weights')

    result = model.fitTo(
        data,
        ROOT.RooFit.AsymptoticError(True),
        # ROOT.RooFit.Minimizer("Minuit2", "minimize"),
        ROOT.RooFit.Save(),
        ROOT.RooFit.PrintLevel(-1),
    )

    frame = x.frame(ROOT.RooFit.Title(f"{name} {fit_type} {sig.name}"))

    return_vals = {
        'nsig': nsig.getVal(),
        'nsig_err': nsig.getError(),
        'nbkg': nbkg.getVal(),
        'mean': sig.mean.getVal(),
        'sigma': sig.sigma.getVal(),
    }

    if fit_type == 'binned':
        return_vals['chi2'] = frame.chiSquare()
    elif fit_type == 'unbinned':
        return_vals['nll'] = result.minNll()

    if plot:
        data.plotOn(frame)
        model.plotOn(frame)

        c = ROOT.TCanvas()
        frame.Draw()

        # Put fit parameters in legend
        legend = ROOT.TLegend(0.1, 0.4, 0.4, 0.9)
        sig.legend(legend)

        if fit_type == 'binned':
            legend.AddEntry(0, f'chi2: {frame.chiSquare():.2f}', '')
        elif fit_type == 'unbinned':
            legend.AddEntry(0, f'nll: {result.minNll():.2f}', '')

        legend.Draw()

        fout=f'{plot_dir}/{name.replace(" ", "_")}_{fit_type}_{sig.name}.png'
        return_vals['filename'] = fout
        c.SaveAs(fout)

    return return_vals

