#! /usr/bin/env python

# modified hangman game -- plays SAT words

import random, sys

hangmanPics = ['''
  +---+
  |   |
      |
      |
      |
      |
=========''', '''
  +---+
  |   |
  O   |
      |
      |
      |
=========''', '''
  +---+
  |   |
  O   |
  |   |
      |
      |
=========''', '''
  +---+
  |   |
  O   |
 /|   |
      |
      |
=========''', '''
  +---+
  |   |
  O   |
 /|\  |
      |
      |
=========''', '''
  +---+
  |   |
  O   |
 /|\  |
 /    |
      |
=========''', '''
  +---+
  |   |
  O   |
 /|\  |
 / \  |
      |
=========''']

#print files
words = {}

rf = open('SAT_words.txt', 'r')
for line in rf.readlines():
  line = line.strip().split()
  word, meaning = line[0].upper(), line[1:]
  words[word] = ' '.join(meaning)
rf.close()

#print words

def getRandomWord(wordDictionary):
  # returns a random word from the dictionary
  word = random.choice(wordDictionary.keys()) # random key / word
  meaning = wordDictionary[word] # meaning for the word
  return word, meaning # return word and meaning

def displayBoard(hangmanPics, missedLetters, correctLetters, secretWord):
  print hangmanPics[len(missedLetters)]

  print '\nMissed letters: ',
  for mletter in missedLetters:
    print mletter,
  print

  blanks = '_' * len(secretWord)

  for i in range(len(secretWord)): # replace blanks with correctly guessed letters
    if secretWord[i] in correctLetters:
      blanks = blanks[:i] + secretWord[i] + blanks[i+1:]

  for letter in blanks: # show the secret word with spaces in between each letter
    print letter,
  print '\n'

def getGuess(alreadyGuessed):
  # returns the letter the player entered. This function makes sure the player
  #   entered a single letter, and not something else.
  while True:
    print 'Guess a letter...'
    guess = raw_input()
    guess = guess.upper()
    if len(guess) != 1:
      print 'Please enter a single letter...'
    elif guess in alreadyGuessed:
      print 'You have already guessed that letter. Choose again ...'
    elif guess not in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
      print 'Please enter a LETTER ...'
    else: return guess

def playAgain():
  # returns True if the player wants to play again, False otherwise
  print '\nDo you want to play again? (yes or no)'
  return raw_input().lower().startswith('y')

print 'H A N G M A N'.center(40)
print 'with SAT/GRE words'.center(40)
print 'written using Python by Laxman Pandey'.center(40)
print '07/19/2014'.center(40)
missedLetters = ''
correctLetters = ''
secretWord, secretWordMeaning = getRandomWord(words)
gameIsDone = False

while True:
  displayBoard(hangmanPics, missedLetters, correctLetters, secretWord)

  # Let the player type in a letter
  guess = getGuess(missedLetters + correctLetters)

  if guess in secretWord:
    correctLetters += guess

    # Check if the player has won
    foundAllLetters = True
    for i in range(len(secretWord)):
      if secretWord[i] not in correctLetters:
        foundAllLetters = False
        break
    if foundAllLetters:
      print 'Yes! The secret word is "%s"! You have won!' % secretWord
      print '%s: %s' % (secretWord, secretWordMeaning)
      gameIsDone = True

  else:
    missedLetters += guess

    # check if player has guessed too many times and lost
    if len(missedLetters) == len(hangmanPics)-1:
      displayBoard(hangmanPics, missedLetters, correctLetters, secretWord)
      print 'You have run out of guesses!'
      print 'You had %s missed guesses and %s correct guesses.' %(
                        len(missedLetters), len(correctLetters))
      print 'The word was "%s"' % secretWord
      print '%s: %s' % (secretWord, secretWordMeaning)
      gameIsDone = True

  # ask if the player wants to play again (but only if the game is done)
  if gameIsDone:
    if playAgain():
      missedLetters = ''
      correctLetters = ''
      gameIsDone = False
      secretWord = getRandomWord(words)
    else:
      print 'Existing on your request ...\n'
      break


