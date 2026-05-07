import json
import os
import streamlit as st

class ISBTDataHandler:
    """Handler for local ISBT blood group allele database with automatic matching."""
    
    BASE_URL = "https://blooddatabase.isbtweb.org"
    ALLELES_DB_PATH = "utils/data/isbt_alleles.json"

    def __init__(self):
        """Initialize the ISBT handler with local allele database."""
        self.alleles_db = self._load_local_alleles()

    def _load_local_alleles(self):
        """Load ISBT allele data from local JSON file."""
        try:
            if os.path.exists(self.ALLELES_DB_PATH):
                with open(self.ALLELES_DB_PATH, 'r') as f:
                    return json.load(f)
            else:
                st.warning(f"Local ISBT database not found at {self.ALLELES_DB_PATH}")
                return {}
        except Exception as e:
            st.warning(f"Error loading ISBT database: {e}")
            return {}

    @st.cache_data(ttl=86400)
    def fetch_alleles(_self, system_code):
        """
        Fetch alleles for a blood group system from local database.
        
        Args:
            system_code: System code like 'ABO', 'RHD', 'KELL', 'DUFFY', 'KIDD', 'HE'
        
        Returns:
            List of allele dictionaries or None if not found.
        """
        if not _self.alleles_db:
            return None
        
        system_data = _self.alleles_db.get(system_code)
        if system_data:
            return list(system_data.values())
        return None

    def get_rhd_variants(self):
        """Get RHD alleles from local database."""
        return self.fetch_alleles("RHD") or []
    
    def get_abo_variants(self):
        """Get ABO alleles from local database."""
        return self.fetch_alleles("ABO") or []
    
    def match_variants_to_alleles(self, variants, system_code):
        """
        Match detected variants to known ISBT alleles for a specific system.
        
        Args:
            variants: List of variant strings (e.g., ["c.261del28", "c.156delA"])
            system_code: Blood group system code (e.g., "ABO", "RHD")
        
        Returns:
            List of matching allele dictionaries
        """
        alleles = self.fetch_alleles(system_code)
        if not alleles:
            return []
        
        matches = []
        for allele in alleles:
            allele_variants = allele.get('variants', [])
            # Check if any detected variant matches this allele
            if any(var in allele_variants for var in variants):
                matches.append(allele)
        
        return matches
    
    def get_system_alleles(self, system_code):
        """Get all alleles for a given system."""
        return self.fetch_alleles(system_code) or []
    
    def suggest_blood_group_from_variants(self, variants_list, system_code):
        """
        Suggest blood group phenotype from detected variants.
        
        Args:
            variants_list: List of detected variant strings
            system_code: Blood group system code
        
        Returns:
            Dictionary with suggested phenotype and matching alleles
        """
        matched_alleles = self.match_variants_to_alleles(variants_list, system_code)
        
        if not matched_alleles:
            return {
                'system': system_code,
                'status': 'no_match',
                'message': 'No matching alleles found in database',
                'matches': []
            }
        
        # Compile phenotype suggestions
        phenotypes = set()
        for allele in matched_alleles:
            phenotypes.add(allele.get('phenotype', 'Unknown'))
        
        return {
            'system': system_code,
            'status': 'matched',
            'phenotypes': list(phenotypes),
            'matches': matched_alleles,
            'message': f"Suggested phenotype: {', '.join(phenotypes)}"
        }
