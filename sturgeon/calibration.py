import numpy as np

class HistogramCalibration():

    def __init__(
        self,
        num_classes = 91,
        num_bins = 10,
    ):

        self.num_classes = num_classes
        self.num_bins = num_bins
        self.conf = np.zeros((self.num_classes, self.num_bins), dtype=float)

        self._calculate_bins()

    def fit(self, probs: np.ndarray, labels: np.ndarray):
        """
        Fit the calibration model, finding optimal confidences for all the bins.
        
        Params:
            probs: probabilities of data [samples, classes]
            labels: numeric labels of the classes
        """
        for label_num in range(self.num_classes):

            subprobs = probs[:, label_num]
            onevs_label = np.zeros(labels.shape, )
            onevs_label[labels == label_num] = 1
        
            for conf_idx, conf_thresh in enumerate(self.upper_bounds):
                acc = onevs_label[(subprobs > conf_thresh-self.bin_size) & (subprobs <= conf_thresh)]
                if len(acc) > 0:
                    acc = np.mean(acc)
                else:
                    acc = conf_thresh
                self.conf[label_num, conf_idx] = acc

    def _calculate_bins(self):
        self.bin_size = 1./self.num_bins
        self.upper_bounds = np.arange(self.bin_size, 1+self.bin_size, self.bin_size)

    def calibrate(self, probs):
        """
        Calibrate a new set of probabilities
        
        Params:
            probs: probabilities of data [classes]
        """
        for i, prob in enumerate(probs):
            idx = np.searchsorted(self.upper_bounds, prob)
            probs[i] = self.conf[i, idx]

        return probs

    def calibrate_batch(self, probs):
        for i in range(probs.shape[0]):
            for j, prob in enumerate(probs[i, :]):
                idx = np.searchsorted(self.upper_bounds, prob)
                probs[i, j] = self.conf[j, idx]

        return probs

    def load(self, calibration_file):
        """
        Args:
            calibration_file (str)
        """

        self.conf = np.load(calibration_file)
        self.num_classes = self.conf.shape[0]
        self.num_bins = self.conf.shape[1]
        self._calculate_bins()

    def load_matrix(self, calibration_matrix):
        """
        Args:
            calibration_file (str)
        """

        self.conf = calibration_matrix
        self.num_classes = self.conf.shape[0]
        self.num_bins = self.conf.shape[1]
        self._calculate_bins()

    def save(self, output_file):

        np.save(output_file, self.conf)
