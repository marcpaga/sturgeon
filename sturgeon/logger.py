import sys
from typing import Optional
import logging

# Logging formatter supporting colorized output
class LogFormatter(logging.Formatter):

	COLOR_CODES = {
		logging.CRITICAL: "\033[1;35m", # bright/bold magenta
		logging.ERROR:    "\033[1;31m", # bright/bold red
		logging.WARNING:  "\033[1;33m", # bright/bold yellow
		logging.INFO:     "\033[0;37m", # white / light gray
		logging.DEBUG:    "\033[1;30m"  # bright/bold black / dark gray
	}

	RESET_CODE = "\033[0m"

	def __init__(self, color, *args, **kwargs):
		super(LogFormatter, self).__init__(*args, **kwargs)
		self.color = color

	def format(self, record, *args, **kwargs):
		if (self.color == True and record.levelno in self.COLOR_CODES):
			record.color_on  = self.COLOR_CODES[record.levelno]
			record.color_off = self.RESET_CODE
		else:
			record.color_on  = ""
			record.color_off = ""
		return super(LogFormatter, self).format(record, *args, **kwargs)

# Setup logging
def setup_logging(
    logfile_file: str, 
    logfile_log_level: Optional[str] = 'debug', 
    console_log_output: Optional[str] = 'stdout',
    console_log_level: Optional[str] = 'info', 
    console_log_color: Optional[bool] = True, 
    log_line_template: Optional[str] = '%(asctime)s - %(name)-12s - %(levelname)-8s - %(message)s',
):


	logger = logging.getLogger()
	logger.setLevel(logging.DEBUG)

	# Create console handler
	console_log_output = console_log_output.lower()
	if (console_log_output == "stdout"):
		console_log_output = sys.stdout
	elif (console_log_output == "stderr"):
		console_log_output = sys.stderr
	else:
		print("Failed to set console output: invalid output: '{}', should be either 'stdout' or 'stderr'".format(console_log_output))
		return False
	console_handler = logging.StreamHandler(console_log_output)

	# Set console log level
	try:
        # only accepts uppercase level names
		console_handler.setLevel(console_log_level.upper()) 
	except:
		print("Failed to set console log level: invalid level: {}".format(console_log_level))
		return False

	# Create and set formatter, add console handler to logger
	console_formatter = LogFormatter(
        fmt = log_line_template, 
        color = console_log_color,
    )
	console_handler.setFormatter(console_formatter)
	logger.addHandler(console_handler)

	# Create log file handler
	try:
		logfile_handler = logging.FileHandler(logfile_file)
	except Exception as exception:
		print("Failed to set up log file: {}".format(exception))
		return False

	# Set log file log level
	try:
        # only accepts uppercase level names
		logfile_handler.setLevel(logfile_log_level.upper()) 
	except:
		print("Failed to set log file log level: invalid level: {}".format(logfile_log_level))
		return False

	# Create and set formatter, add log file handler to logger
	logfile_formatter = LogFormatter(
        fmt = log_line_template, 
        color = False,
    )
	logfile_handler.setFormatter(logfile_formatter)
	logger.addHandler(logfile_handler)

	# Success
	return True