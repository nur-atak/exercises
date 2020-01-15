-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='ONLY_FULL_GROUP_BY,STRICT_TRANS_TABLES,NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION';

-- -----------------------------------------------------
-- Schema exercise4
-- -----------------------------------------------------

-- -----------------------------------------------------
-- Schema exercise4
-- -----------------------------------------------------
CREATE SCHEMA IF NOT EXISTS `exercise4` ;
USE `exercise4` ;

-- -----------------------------------------------------
-- Table `exercise4`.`Doctor`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `exercise4`.`Doctor` (
  `idDoctor` INT NOT NULL,
  `Firstname` VARCHAR(45) NULL,
  `Lastname` VARCHAR(45) NULL,
  `Date_of_birth` DATE NULL,
  `Address` BLOB NULL,
  `Phonenumber` VARCHAR(20) NULL,
  PRIMARY KEY (`idDoctor`),
  UNIQUE INDEX `idDoctor_UNIQUE` (`idDoctor` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `exercise4`.`Medical`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `exercise4`.`Medical` (
  `idMedical` INT NOT NULL,
  `Overtime_rate` VARCHAR(45) NULL,
  `idDoctor` INT NULL,
  PRIMARY KEY (`idMedical`),
  INDEX `idDoctor_idx` (`idDoctor` ASC) VISIBLE,
  UNIQUE INDEX `idMedical_UNIQUE` (`idMedical` ASC) VISIBLE,
  CONSTRAINT `idDoctor1`
    FOREIGN KEY (`idDoctor`)
    REFERENCES `exercise4`.`Doctor` (`idDoctor`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `exercise4`.`Specialist`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `exercise4`.`Specialist` (
  `idSpecialist` INT NOT NULL,
  `Field_area` VARCHAR(45) NULL,
  `idDoctor` INT NULL,
  PRIMARY KEY (`idSpecialist`),
  UNIQUE INDEX `idSpecialist_UNIQUE` (`idSpecialist` ASC) VISIBLE,
  INDEX `idDoctor_idx` (`idDoctor` ASC) VISIBLE,
  CONSTRAINT `idDoctor2`
    FOREIGN KEY (`idDoctor`)
    REFERENCES `exercise4`.`Doctor` (`idDoctor`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `exercise4`.`Patient`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `exercise4`.`Patient` (
  `idPatient` INT NOT NULL,
  `Name` VARCHAR(45) NULL,
  `Address` BLOB NULL,
  `Phone_number` VARCHAR(20) NULL,
  `Date_of_birth` DATE NULL,
  PRIMARY KEY (`idPatient`),
  UNIQUE INDEX `idPatient_UNIQUE` (`idPatient` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `exercise4`.`Appointment`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `exercise4`.`Appointment` (
  `idAppointment` INT NOT NULL,
  `Date` DATE NULL,
  `Time` VARCHAR(45) NULL,
  `Patient` INT NULL,
  `Doctor` INT NULL,
  PRIMARY KEY (`idAppointment`),
  INDEX `Patient_idx` (`Patient` ASC) VISIBLE,
  INDEX `Doctor_idx` (`Doctor` ASC) VISIBLE,
  UNIQUE INDEX `idAppointment_UNIQUE` (`idAppointment` ASC) VISIBLE,
  CONSTRAINT `Patient4`
    FOREIGN KEY (`Patient`)
    REFERENCES `exercise4`.`Patient` (`idPatient`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `Doctor4`
    FOREIGN KEY (`Doctor`)
    REFERENCES `exercise4`.`Doctor` (`idDoctor`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `exercise4`.`Payment`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `exercise4`.`Payment` (
  `idPayment` INT NOT NULL,
  `Method` VARCHAR(45) NULL,
  `Details` BLOB NULL,
  PRIMARY KEY (`idPayment`),
  UNIQUE INDEX `idPayment_UNIQUE` (`idPayment` ASC) VISIBLE,
  UNIQUE INDEX `Method_UNIQUE` (`Method` ASC) VISIBLE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `exercise4`.`Bill`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `exercise4`.`Bill` (
  `idBill` INT NOT NULL,
  `Patient` INT NULL,
  `Doctor` INT NULL,
  `Total` INT NULL,
  PRIMARY KEY (`idBill`),
  UNIQUE INDEX `idBill_UNIQUE` (`idBill` ASC) VISIBLE,
  INDEX `Patient_idx` (`Patient` ASC) VISIBLE,
  INDEX `Doctor_idx` (`Doctor` ASC) VISIBLE,
  CONSTRAINT `Patient`
    FOREIGN KEY (`Patient`)
    REFERENCES `exercise4`.`Patient` (`idPatient`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `Doctor`
    FOREIGN KEY (`Doctor`)
    REFERENCES `exercise4`.`Doctor` (`idDoctor`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `exercise4`.`Transaction`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `exercise4`.`Transaction` (
  `idTransaction` INT UNSIGNED NOT NULL,
  `Payment` INT NULL,
  `Bill` INT NULL,
  PRIMARY KEY (`idTransaction`),
  INDEX `Payment_idx` (`Payment` ASC) VISIBLE,
  INDEX `Bill_idx` (`Bill` ASC) VISIBLE,
  CONSTRAINT `Payment3`
    FOREIGN KEY (`Payment`)
    REFERENCES `exercise4`.`Payment` (`idPayment`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `Bill3`
    FOREIGN KEY (`Bill`)
    REFERENCES `exercise4`.`Bill` (`idBill`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;
