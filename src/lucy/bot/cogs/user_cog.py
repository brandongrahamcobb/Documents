''' user_cog.py
    Copyright (C) 2024 github.com/brandongrahamcobb

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

from bot.main import Lucy
from discord import app_commands
from discord.ext import commands
from rdkit import Chem
from rdkit.Chem import Crippen
from typing import List

import bot.utils.helpers as lucy

import discord
import io
import os
import pubchempy as pcp
import rdkit
import requests
import shlex
import traceback

class UserCog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot
        self.config = bot.config
        self.stacks = {}
        self.diaminobenzidine = lucy.get_mol('diaminobenzidine')

    def get_user_stack(self, user_id: int):
        return self.stacks.get(user_id, [])

    def set_user_stack(self, user_id: int, molecule_names: List[str]):
        mol_objects = [lucy.get_mol(name) for name in molecule_names]
        self.stacks[user_id] = mol_objects

    @commands.command(name='stack', description='Chose your guide.')
    async def analog(self, ctx: commands.Context, *args):
        new_stack = list(args)
        if args:
            self.set_user_stack(user_id=ctx.author.id, molecule_names=args)
        original = max(self.get_user_stack(ctx.author.id), key=lambda mol: lucy.get_proximity(mol, self.diaminobenzidine))
        while True:
            await ctx.send(f'{ctx.author.name}\'s Stack')
            await ctx.send([lucy.get_molecule_name(mol) for mol in self.get_user_stack(ctx.author.id)])
            await ctx.send('Which molecule would you like to interview for your stack?')
            response = await self.bot.wait_for(
                'message',
                timeout=20000.0,
                check=lambda message: message.author == ctx.author and message.channel == ctx.channel
            )
            try:
                new = lucy.get_mol(response.content)
                original_proximity = lucy.get_proximity(self.diaminobenzidine, original)
                new_proximity = lucy.get_proximity(self.diaminobenzidine, new)
                if new_proximity < original_proximity:
                    new_stack.append(response.content)
                    embed = lucy.get_sds(response.content)
                    await ctx.send(embed=embed)
                    self.set_user_stack(user_id=ctx.author.id, molecule_names=new_stack)
                elif new_proximity > original_proximity:
                    await ctx.send(f'DANGER: {response.content}')
                await ctx.send(f'{(100 * new_proximity):.3f}% {response.content} {(100 * original_proximity):.3f}%')
            except Exception as e:
                await ctx.send(e)
                break

    @commands.hybrid_command(name='compare', description='Draw a molecule or compare molecules by their names.')
    async def compare(self, ctx: commands.Context, *, molecules: str) -> None:
        if ctx.interaction:
            await ctx.interaction.response.defer(ephemeral=True)
        args = shlex.split(molecules)
        mol1 = lucy.get_mol(args[0])
        mol2 = lucy.get_mol(args[1])
        image = lucy.view_difference(mol1, mol2)
        await ctx.send(file=image)

    @commands.hybrid_command(name='draw', description='Draw a molecule or compare molecules by their names.')
    async def draw(self, ctx: commands.Context, *, molecules: str) -> None:
        try:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            args = shlex.split(molecules)
            if len(args) == 1:
                mol = lucy.get_mol(args[0])
                if mol is None:
                    embed = 'Invalid molecule name or structure.'
                    await ctx.send(embed=embed)
                    return
                image = lucy.draw_watermarked_molecule(mol)
                await ctx.send(file=discord.File(image, f'{args[0]}.png'))
                return
            pairs = lucy.unique_pairs(args)
            if not pairs:
                embed = 'No valid pairs found.'
                await ctx.send(embed=embed)
                return
            for pair in pairs:
                mol = lucy.get_mol(pair[0])
                refmol = lucy.get_mol(pair[1])
                if mol is None or refmol is None:
                    embed = f'One or both of the molecules {pair[0]} or {pair[1]} are invalid.'
                    await ctx.send(embed=embed)
                    continue
                #image = lucy.draw_watermarked_molecule(mol)
                #proximity = lucy.get_proximity(default=mol, input=refmol)
                #embed = f'The Tanimoto Similarity of {lucy.get_molecule_name(mol)} and {lucy.get_molecule_name(refmol)} is {proximity:.2f}'
                #await interaction.followup.send(embed=embed)
                fingerprints = [
                   lucy.draw_fingerprint([mol, refmol]),
                    lucy.draw_fingerprint([refmol, mol])
                ]
                combined_image = lucy.combine(fingerprints[0], pair[1], fingerprints[1], pair[0])
                await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))
        except Exception as e:
            await ctx.send(f'An error occurred: {str(e)}')

    @commands.hybrid_command(name='emoji')
    async def emoji(self, ctx: commands.Context, *, argument: str):
        if ctx.interaction:
             await ctx.interaction.response.defer(ephemeral=True)
        await ctx.send(lucy.get_emoji(argument))

    @commands.hybrid_command(name='get', hidden=True)
    async def get(self, ctx: commands.Context, *, argument: str):
#        if ctx.guild.id != 1131418877214068858:
 #           return
        if ctx.interaction:
             await ctx.interaction.response.defer(ephemeral=True)
        try:
            file = lucy.stable_cascade(argument)
            if isinstance(file, discord.File):
                await ctx.send(file=file)
            else:
                await ctx.send(f"Error generating image: {file}")
        except Exception as e:
            print(f"Error in on_message: {e}")
            await ctx.send(f"An unexpected error occurred: {e}")

    @commands.hybrid_command(name='sim')
    async def sim(self, ctx: commands.Context, *, molecules: str):
        args = shlex.split(molecules)
        similarity = lucy.get_proximity(lucy.get_mol(args[0]), lucy.get_mol(args[1]))
        await ctx.send(similarity)

    @commands.hybrid_command(name='logp')
    async def logp(self, ctx: commands.Context, *, molecules: str):
        args = shlex.split(molecules)
        for arg in args:
            compounds = pcp.get_compounds(arg, 'name')
            compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
            mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
            log_p = Crippen.MolLogP(mol)
            await ctx.send(f'Your octanol:water coefficient is: {log_p}')

    @commands.hybrid_command(name='msds')
    async def msds(self, ctx: commands.Context, *, argument: str):
       embed = lucy.get_sds(argument)
       await ctx.send(embed=embed)

    @commands.hybrid_command(name='search')
    async def search(self, ctx: commands.Context, *, molecules: str):
        try:
            if ctx.interaction:
                 await ctx.interaction.response.defer(ephemeral=True)
            args = shlex.split(molecules)
            try:
                for arg in args:
                    watermarked_image = lucy.gsrs(molecules)
                    with io.BytesIO() as image_binary:
                        watermarked_image.save(image_binary, format='PNG')
                        image_binary.seek(0)
                        await ctx.send(file=discord.File(fp=image_binary, filename='watermarked_image.png'))
            except:
                await ctx.send(traceback.format_exc())
        except:
            if arg is None:
                return
            await ctx.send(f'{arg} is an unknown molecule.')

    @commands.hybrid_command(name='snort')
    async def snort(self, ctx: commands.Context, molecules: str, quantity: int):
        fingerprints = []
        names = []
        if not molecules or molecules == 'cocaine' or 'cocaethylene' or 'coke' and quantity > 100:
            if ctx.interaction:
                await ctx.interaction.response.defer(ephemeral=True)
            args = shlex.split(molecules)
            molecule = lucy.get_mol(args[0])
            for mol in range(quantity):
                names.append(args[0])
                fingerprints.append(lucy.draw_fingerprint([molecule, molecule]))
            combined_image = lucy.combine(fingerprints, names)
            await ctx.send(file=discord.File(combined_image, f'molecule_comparison.png'))

    @commands.hybrid_command(name='smiles')
    async def smiles(self, ctx: commands.Context, *, molecules: str):
        try:
            if ctx.interaction:
                 await ctx.interaction.response.defer(ephemeral=True)
            args = shlex.split(molecules)
            try:
                for arg in args:
                    compounds = pcp.get_compounds(arg, 'name')
                    compound = compounds[0]
                    isomeric_smiles = compound.isomeric_smiles
                    await ctx.send(f'The isomeric SMILES for {arg} is: {isomeric_smiles}')
            except:
                await ctx.send(traceback.format_exc())
        except:
            if arg is None:
                return
            await ctx.send(f'{arg} is an unknown molecule.')

async def setup(bot: commands.Bot):
    await bot.add_cog(UserCog(bot))
