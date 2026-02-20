

"""Configuration handler for TranslonScorer."""
import os
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union, Any

@dataclass
class Config:
    """Configuration class to handle and validate TranslonScorer parameters."""
    
    # Input files
    bam: Optional[str] = None
    bam_collapsed: bool = False
    bam_count_from: Optional[str] = None  # 'name' or 'tag'
    bam_count_pattern: Optional[str] = None
    bam_count_tag: Optional[str] = None
    chromsizes: Optional[str] = None
    sequence: str = field(default="")
    annotation: str = field(default="")
    
    # BigWig options
    bigwig: Optional[str] = None
    forward_bigwig: Optional[str] = None
    reverse_bigwig: Optional[str] = None

    # Zarr unique read matrix
    zarr_root: Optional[str] = None
    read_index_parquet: Optional[str] = None
    samples: Optional[List[str]] = None
    
    # Analysis options
    stranded: bool = False
    offsets: Optional[str] = None
    start_codons: List[str] = field(default_factory=lambda: ["ATG"])
    stop_codons: List[str] = field(default_factory=lambda: ["TAA", "TAG", "TGA"])
    min_length: int = 0
    max_length: int = 1000000
    sru_range: int = 15
    scoring_method: str = "modern"
    plot_range: int = 30
    
    # Output options
    output: str = field(default="")
    log_file: Optional[str] = None
    log_level: str = "INFO"
    
    # Runtime state
    bigwig_paths: Union[str, Dict[str, str], None] = None
    
    @classmethod
    def from_click_args(cls, **kwargs) -> 'Config':
        """Create a Config instance from Click command arguments."""
        # Process any string lists that come as comma-separated values
        if 'start_codons' in kwargs and isinstance(kwargs['start_codons'], str):
            kwargs['start_codons'] = kwargs['start_codons'].split(',')
        
        if 'stop_codons' in kwargs and isinstance(kwargs['stop_codons'], str):
            kwargs['stop_codons'] = kwargs['stop_codons'].split(',')
        
        # Handle naming differences between CLI and config
        if 'bam_path' in kwargs:
            kwargs['bam'] = kwargs.pop('bam_path')
        
        if 'bigwig_path' in kwargs:
            kwargs['bigwig'] = kwargs.pop('bigwig_path')
            
        if 'outfile' in kwargs:
            kwargs['output'] = kwargs.pop('outfile')
            
        # Filter out None values for optional parameters to use defaults
        filtered_kwargs = {k: v for k, v in kwargs.items() if v is not None or k in ['stranded']}
        
        return cls(**filtered_kwargs)
    
    def validate(self) -> None:
        """Validate the configuration and set derived values."""
        # Check required fields
        if not self.sequence:
            raise ValueError("Sequence file is required")
        
        if not self.annotation:
            raise ValueError("Annotation file is required")
        
        if not self.output:
            raise ValueError("Output path is required")
            
        # Validate input files exist
        self._validate_file_exists(self.sequence, "Sequence")
        self._validate_file_exists(self.annotation, "Annotation")
        
        if self.bam:
            self._validate_file_exists(self.bam, "BAM")
            if not self.chromsizes and not (self.bigwig or (self.forward_bigwig and self.reverse_bigwig)):
                raise ValueError("Chromosome sizes file is required when processing BAM files")
                
        if self.chromsizes:
            self._validate_file_exists(self.chromsizes, "Chromosome sizes")
            
        if self.offsets:
            self._validate_file_exists(self.offsets, "Offsets")
            
        # Validate BigWig files
        if self.bigwig:
            self._validate_file_exists(self.bigwig, "BigWig")
        
        if self.forward_bigwig:
            self._validate_file_exists(self.forward_bigwig, "Forward BigWig")
            if not self.reverse_bigwig:
                raise ValueError("Reverse BigWig must be provided with Forward BigWig")
                
        if self.reverse_bigwig:
            self._validate_file_exists(self.reverse_bigwig, "Reverse BigWig")
            if not self.forward_bigwig:
                raise ValueError("Forward BigWig must be provided with Reverse BigWig")
        
        # Validate input combinations
        if not any([self.bam, self.bigwig, (self.forward_bigwig and self.reverse_bigwig), self.zarr_root]):
            raise ValueError(
                "Provide one of: BAM, BigWig, forward+reverse BigWigs, or Zarr root"
            )
            
        # Set up bigwig_paths based on inputs
        self._setup_bigwig_paths()
        
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(self.output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    def _validate_file_exists(self, filepath: str, description: str) -> None:
        """Check if a file exists and raise an error if it doesn't."""
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"{description} file not found: {filepath}")
    
    def _setup_bigwig_paths(self) -> None:
        """Set up bigwig_paths based on the provided inputs."""
        # Handle strand-specific bigwig inputs
        if self.forward_bigwig and self.reverse_bigwig:
            self.bigwig_paths = {
                'forward': self.forward_bigwig, 
                'reverse': self.reverse_bigwig
            }
            self.stranded = True
        elif self.bigwig:
            self.bigwig_paths = self.bigwig
        else:
            self.bigwig_paths = None
            
    def setup_logging(self) -> None:
        """Configure logging based on the provided log level and file."""
        log_level = getattr(logging, self.log_level)
        
        # Create basic configuration
        logging_config = {
            'level': log_level,
            'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            'datefmt': '%Y-%m-%d %H:%M:%S'
        }
        
        # Add file handler if log_file is specified
        if self.log_file:
            log_dir = os.path.dirname(self.log_file)
            if log_dir and not os.path.exists(log_dir):
                os.makedirs(log_dir)
            logging_config['filename'] = self.log_file
        
        # Apply configuration
        logging.basicConfig(**logging_config)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the configuration to a dictionary for serialization."""
        # Get all fields as a dictionary
        config_dict = {
            key: getattr(self, key) 
            for key in self.__dataclass_fields__.keys()
        }
        
        # Remove internal/runtime fields
        if 'bigwig_paths' in config_dict:
            config_dict.pop('bigwig_paths')
        
        return config_dict
