{- 
  ribosome-2.hs - A high-level model of the action of a ribosome in creating a polypeptide from an
                  mRNA strand.  In reality this is much more detailed, I'll add further detail in
                  subsequent revisions.

                  The prior file ribosome.hs used the first part of the finite state machine blog
                  to implement the synthesis process. This revision uses the second part, which
                  separates the protocol from the program. This is the first step in generalizing
                  the idea of synthesizing something by traversing a precursor product.

  References: 

    Wiki entry on ribosomes
      https://en.wikipedia.org/wiki/Ribosome

    Video showing the polypeptide construction process, from DNA Learning Center 
      https://www.youtube.com/watch?v=TfYf_rPWUdY

    Finite state machines in Haskell, described by Oskar Wickström
      https://wickstrom.tech/finite-state-machines/2017/11/10/finite-state-machines-part-1-modeling-with-haskell.html
      https://wickstrom.tech/finite-state-machines/2017/11/19/finite-state-machines-part-2.html

    Map of codons to amino acids
      https://www.genome.jp/kegg/catalog/codes1.html
-}

{-# LANGUAGE GADTs                      #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TypeFamilies               #-}

import Control.Monad.IO.Class

{- Set up our types -}

-- ACGU
type Nucleotide  = Char
type Codon       = (Nucleotide, Nucleotide, Nucleotide)

-- A new mRNA strand has introns interspersed with exons, but for now we'll assume the introns have
-- been spliced out, though this is a precursor process done by a spliceosome, we'll model it later
type MRNA        = [Nucleotide]

-- A polypeptide is a sequence of linked amino acids emitted from the ribosome and folded later
-- into a protein. Building this from an mRNA strand is the purpose of the code in this file
type Polypeptide = [AminoAcid]


{- Below are roughly the states a polypeptide assembly can be in. They are "phantom types" that
   carry no data, they serve only as markers for the state machine (see blog linked above) -}

-- An mRNA strand is floating around the cytosol, ready to be processed by a ribosome
data MRNAOnly

-- The ribosome has been assembled from its constituent subunits (40S, 60S) and can process codons
data RibosomeLoaded

-- Amino acid assembly is complete, a polypeptide is born
data Assembled

-- Our concrete state type, holding the accumulated data as function args
data SynthesisState s where
  MRNAOnly                :: MRNA                -> SynthesisState MRNAOnly
  RibosomeLoaded          :: MRNA -> Polypeptide -> SynthesisState RibosomeLoaded
  Assembled               :: Polypeptide         -> SynthesisState Assembled


-- Encode the events of the synthesis pathway as a type class. This is a bit specialized to protein
-- synthesis, but we'll generalize it in subsequent revisions
class ProteinSynthesis m where
  type State m :: * -> *

  initial         :: MRNA                      -> m (State m MRNAOnly)
  recruitRibosome :: (State m MRNAOnly)        -> m (State m RibosomeLoaded)
  processCodon    :: (State m RibosomeLoaded)  -> m (State m RibosomeLoaded)
  isFinished      :: (State m RibosomeLoaded)  -> m Bool
  detach          :: (State m RibosomeLoaded)  -> m (State m Assembled)
  end             :: (State m Assembled)       -> m Polypeptide


-- The type we'll instantiate our ProteinSynthesis class on
newtype ProteinSynthesisT m a = ProteinSynthesisT { runProteinSynthesisT :: m a }
                                deriving (Functor, Monad, Applicative, MonadIO)

instance (MonadIO m) => ProteinSynthesis (ProteinSynthesisT m) where
  type State (ProteinSynthesisT m) = SynthesisState

  initial mRNA                    = return $ MRNAOnly mRNA
  recruitRibosome (MRNAOnly mRNA) = return $ RibosomeLoaded mRNA []

  processCodon (RibosomeLoaded mRNA poly)
                                  = return $ RibosomeLoaded mRNA' poly'
                                             where mRNA' = drop 3 mRNA
                                                   codon = nextCodon mRNA
                                                   poly' = toAmino codon : poly

  isFinished (RibosomeLoaded mRNA _) = return $ empty || stop
                                                where empty = null mRNA
                                                      stop  = isStopCodon $ nextCodon mRNA

  detach (RibosomeLoaded _ poly)     = return $ Assembled (reverse poly)
  end (Assembled poly)               = return poly
      

-- The program itself, it drives the state machine
synthesize :: (ProteinSynthesis m, MonadIO m) => MRNA -> m Polypeptide
synthesize mRNA = initial mRNA
                  >>= recruitRibosome
                  >>= processCodons
                  >>= end

-- Recursively process the next codon until done
processCodons :: (ProteinSynthesis m, MonadIO m) => State m RibosomeLoaded -> m (State m Assembled)
processCodons s = do finished <- isFinished s
                     if finished
                       then detach s
                       else processCodon s 
                            >>= processCodons

assemble :: MRNA -> IO ()
assemble mRNA = do
  polypeptide <- runProteinSynthesisT $ synthesize mRNA
  putStrLn $ "Assembled polypeptide: " ++ show polypeptide


{- Test run using the mRNA strand listed in the wiki picture at
   https://en.wikipedia.org/wiki/Ribosome#/media/File:Peptide_syn.png

    *Main> assemble "UGGAAAGAUUUCUAGUUCUUC"  -- Note the stop codon UAG I've added
    Assembled polypeptide: [Trp,Lys,Asp,Phe]
-}

-- Amino acid data, this will be factored into its own file eventually
data AminoAcid = Ala | Arg | Asn | Asp | Cys
               | Glu | Gln | Gly | His | Ile
               | Leu | Lys | Met | Phe | Pro
               | Ser | Thr | Trp | Tyr | Val
               deriving (Show)

nextCodon :: MRNA -> Codon
nextCodon (n1:n2:n3:_) = (n1, n2, n3)

isStopCodon :: Codon -> Bool
isStopCodon (n1,n2,n3) = [n1,n2,n3] `elem` ["UAG", "UAA", "UGA"]

-- Map a triplet of mRNA nucleotides (a codon) to its corresponding amino acid
toAmino :: Codon -> AminoAcid
toAmino (n1,n2,n3) = case [n1,n2,n3] of 
  { 
    "AAA" -> Lys; "AAG" -> Lys;
    "AAC" -> Asn; "AAU" -> Asn;
    "ACA" -> Thr; "ACC" -> Thr; "ACG" -> Thr; "ACU" -> Thr;
    "AGA" -> Arg; "AGG" -> Arg; "CGA" -> Arg; "CGC" -> Arg; "CGG" -> Arg; "CGU" -> Arg;
    "AGC" -> Ser; "AGU" -> Ser; "UCA" -> Ser; "UCC" -> Ser; "UCG" -> Ser; "UCU" -> Ser;
    "AUA" -> Ile; "AUC" -> Ile; "AUU" -> Ile;
    "AUG" -> Met;
    "CAA" -> Gln; "CAG" -> Gln;
    "CAC" -> His; "CAU" -> His;
    "CCA" -> Pro; "CCC" -> Pro; "CCG" -> Pro; "CCU" -> Pro;
    "CUA" -> Leu; "CUC" -> Leu; "CUG" -> Leu; "CUU" -> Leu; "UUA" -> Leu; "UUG" -> Leu;
    "GAA" -> Glu; "GAG" -> Glu;
    "GAC" -> Asp; "GAU" -> Asp;
    "GCA" -> Ala; "GCC" -> Ala; "GCG" -> Ala; "GCU" -> Ala;
    "GGA" -> Gly; "GGC" -> Gly; "GGG" -> Gly; "GGU" -> Gly;
    "GUA" -> Val; "GUC" -> Val; "GUG" -> Val; "GUU" -> Val;
    "UAC" -> Tyr; "UAU" -> Tyr;
    "UGC" -> Cys; "UGU" -> Cys;
    "UGG" -> Trp;
    "UUC" -> Phe; "UUU" -> Phe;
    _     -> error "Unknown codon";
  }

