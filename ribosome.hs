{- 
  ribosome.hs - A high-level model of the action of a ribosome in creating a
                polypeptide from an mRNA strand.  In reality this is much more
                detailed, I'll add further detail in subsequent revisions

  References: 

    Wiki entry on ribosomes
    https://en.wikipedia.org/wiki/Ribosome

    Video showing the polypeptide construction process, from DNA Learning Center 
    https://www.youtube.com/watch?v=TfYf_rPWUdY

    Finite state machines in Haskell, described by Oskar Wickström
    https://wickstrom.tech/finite-state-machines/2017/11/10/finite-state-machines-part-1-modeling-with-haskell.html

    Map of codons to aminos acid
    https://www.genome.jp/kegg/catalog/codes1.html
-}

import Control.Monad (foldM)
import Text.Printf   (printf)

-- ACGU
type Nucleotide  = Char

-- Messenger RNA is mapped from the DNA strand during transcription. A new mRNA
-- strand has introns interspersed with exons, but for now we'll assume the
-- introns have been spliced out, though this is a precursor process done by
-- spliceosomes, which we can model later
type MRNA        = [Nucleotide]

-- A polypeptide is a sequence of linked amino acids emitted from the ribosome
-- and folded later into a protein. Building this from a starting mRNA strand
-- is the purpose of the code in file
type Polypeptide = [AminoAcid]

-- Roughly the states a polypeptide assembly can be in and the "data" required
-- to describe the state
data ProteinSynthesisState = 

  -- A mRNA strand is floating around the cytosol, ready to be processed by
  -- a ribosome
  MRNAOnly MRNA

  -- The ribosome has been assembled from its constituent subunits (40S,60S)
  -- and is loaded onto the start codon of the mRNA strand, ready to go. We
  -- don't yet model the start/stop codons
  | RibosomeLoaded MRNA

  -- The main work of the ribosome, collecting appropriate tRNA molecules
  -- loaded with the amino acid required next
  | RibosomeWithPolypeptide MRNA Polypeptide

  -- Amino acid assembly is complete, a polypeptide is born
  | Assembled Polypeptide
  deriving (Show)

-- Events in the polypeptide assembly workflow
data ProteinSynthesisEvent
  = RecruitRibosome
  | ProcessCodon
  | Detach
  deriving (Show)

-- Finite state machine (see the blog entry linked in the references above)
type FSM s e = s -> e -> IO s

-- The various transitions between states given a state and event
synthesize :: FSM ProteinSynthesisState ProteinSynthesisEvent

synthesize (MRNAOnly mRNA) RecruitRibosome = 
  return $ RibosomeLoaded mRNA

synthesize (RibosomeLoaded mRNA) ProcessCodon =
  return $ RibosomeWithPolypeptide mRNA' polypeptide'
           where polypeptide' = [amino]
                 amino        = toAmino $ take 3 mRNA
                 mRNA'        = drop 3 mRNA

synthesize (RibosomeWithPolypeptide mRNA polypeptide) ProcessCodon =
  return $ RibosomeWithPolypeptide mRNA' polypeptide'
           where polypeptide' = polypeptide ++ [amino]
                 amino        = toAmino $ take 3 mRNA
                 mRNA'        = drop 3 mRNA

synthesize (RibosomeWithPolypeptide mRNA polypeptide) Detach =
  if null mRNA 
    then return $ Assembled polypeptide
    else error "Detaching early from mRNA strand, assembly incomplete"

-- Run the finite machine by folding over a list of events
runFsm :: Foldable f => FSM s e -> s -> f e -> IO s
runFsm = foldM

-- Borrowed verbatim from the blog linked in the references above, this adds
-- text printing to the state machine to track the evolution of the process
withLogging :: (Show s, Show e) => FSM s e -> FSM s e
withLogging fsm s e = do
  s' <- fsm s e
  printf "- %s x %s -> %s\n" (show s) (show e) (show s')
  return s'

-- Assume a 4-codon mRNA (hence 4 ProcessCodon events)
assemble :: MRNA -> IO ProteinSynthesisState
assemble mRNA = runFsm (withLogging synthesize) start events
  where start  = MRNAOnly mRNA
        events = [ RecruitRibosome,
                   ProcessCodon,
                   ProcessCodon,
                   ProcessCodon,
                   ProcessCodon,
                   Detach ]

{- Test run using the mRNA strand listed in the wiki picture at
   https://en.wikipedia.org/wiki/Ribosome#/media/File:Peptide_syn.png

    *Main> assemble "UGGAAAGAUUUC"
      - MRNAOnly "UGGAAAGAUUUC" x RecruitRibosome -> RibosomeLoaded "UGGAAAGAUUUC"
      - RibosomeLoaded "UGGAAAGAUUUC" x ProcessCodon -> RibosomeWithPolypeptide "AAAGAUUUC" [Trp]
      - RibosomeWithPolypeptide "AAAGAUUUC" [Trp] x ProcessCodon -> RibosomeWithPolypeptide "GAUUUC" [Trp,Lys]
      - RibosomeWithPolypeptide "GAUUUC" [Trp,Lys] x ProcessCodon -> RibosomeWithPolypeptide "UUC" [Trp,Lys,Asp]
      - RibosomeWithPolypeptide "UUC" [Trp,Lys,Asp] x ProcessCodon -> RibosomeWithPolypeptide "" [Trp,Lys,Asp,Phe]
      - RibosomeWithPolypeptide "" [Trp,Lys,Asp,Phe] x Detach -> Assembled [Trp,Lys,Asp,Phe]
      Assembled [Trp,Lys,Asp,Phe]
-}

-- Amino acid data, this will be factored into its own file eventually
data AminoAcid = Ala | Arg | Asn | Asp | Cys
               | Glu | Gln | Gly | His | Ile
               | Leu | Lys | Met | Phe | Pro
               | Ser | Thr | Trp | Tyr | Val
               deriving (Show)

-- Map a tiplet of mRNA nucleotides (a codon) to its corresponding amino acid
toAmino :: [Nucleotide] -> AminoAcid
toAmino ns = case ns of 
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

